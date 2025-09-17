function [Tddot_cmd, tau_roll, tau_pitch, tau_yaw, gimRollAccCmd, gimPitchAccCmd, ghost] = ...
    M_block(u_fw, qc_state, params, varargin)
% M_BLOCK (robust, paper-aligned, seatbelt-toggleable)
%   Maps FW inputs to QC+gimbal commands via DFL/NDI using an internal ghost FW.
%
% Inputs:
%   u_fw     : [thrust; elevator; aileron; rudder]
%   qc_state : [x y z vN vE vD q0 q1 q2 q3 p q r T nu_T eta theta_g eta_dot theta_g_dot]'
%   params   : struct (see below)
%   varargin : optional seatbelt flag (logical or "on"/"off"), default = ON
%
% Outputs:
%   Tddot_cmd, tau_roll, tau_pitch, tau_yaw, gimRollAccCmd, gimPitchAccCmd, ghost
%
% Key params / defaults:
%   dt,g,m_qc,J_qc,tau_max
%   k0_pos..k3_pos, kpsi0,kpsi1
%   dfl.lambda0, dfl.lambda_b3, dfl.lambda_T, dfl.Tref, dfl.T_floor_frac
%   ghost_alpha_acc, ghost_alpha_jerk
%   kT0,kT1, qc_Tddot_max
%   caps.s_max (3x1), caps.psidd_max (scalar)
%   fb.kR(3), fb.kOm(3), fb.theta_act, fb.max_blend, fb.enableSeatbelt (default true)
%
% Diagnostics (ghost.flags):
%   quatReset, nanClamped, illCondSolve, tauSaturated, gainShapeFix, demandClamped

% ---------- unpack QC (columnize)
qc_state = qc_state(:);
xq=qc_state(1); yq=qc_state(2); zq=qc_state(3);
vN=qc_state(4); vE=qc_state(5); vD=qc_state(6);
qBN_qc=qc_state(7:10);
p=qc_state(11); q=qc_state(12); r=qc_state(13);
T=qc_state(14); nu_T=qc_state(15);
eta=qc_state(16); thg=qc_state(17); etad=qc_state(18); thgd=qc_state(19);

% ---------- params & defaults
dt  = getdef(params,'dt',0.001);
g   = getdef(params,'g',9.81);
m   = getdef(params,'m_qc',1.0);
J   = getdef(params,'J_qc',diag([0.02 0.02 0.035]));
tau_lim = getdef(params,'tau_max',[1.2;1.2;1.2]);

% outer-loop gains (can be scalars or 3-vectors)
k0 = getdef(params,'k0_pos',3.0);
k1 = getdef(params,'k1_pos',6.0);
k2 = getdef(params,'k2_pos',4.0);
k3 = getdef(params,'k3_pos',2.0);

kpsi0 = getdef(params,'kpsi0',3.0);
kpsi1 = getdef(params,'kpsi1',1.5);

% regularization
lam0   = getdef_nested(params,'dfl','lambda0',1e-3);
lam_b3 = getdef_nested(params,'dfl','lambda_b3',5e-2);
lam_T  = getdef_nested(params,'dfl','lambda_T',0.8);
Tref   = getdef_nested(params,'dfl','Tref',5);

% ghost smoothing
alpha_a = getdef(params,'ghost_alpha_acc',0.4);
alpha_j = getdef(params,'ghost_alpha_jerk',0.4);

% thrust channel
kT0 = getdef(params,'kT0',80);
kT1 = getdef(params,'kT1',20);
Tdd_lim = getdef(params,'qc_Tddot_max',Inf);

% caps
s_max     = getdef_nested(params,'caps','s_max',[8;8;8]);
psidd_max = getdef_nested(params,'caps','psidd_max',0.50);

% gimbal PD gains
kg0_r=getdef(params,'kg0_roll',80);  kg1_r=getdef(params,'kg1_roll',20);
kg0_p=getdef(params,'kg0_pitch',80); kg1_p=getdef(params,'kg1_pitch',20);

% safety attitude feedback (seatbelt)
kR       = getdef_nested(params,'fb','kR',[0.25;0.25;0.35]);
kOm      = getdef_nested(params,'fb','kOm',[0.03;0.03;0.04]);
thetaAct = getdef_nested(params,'fb','theta_act',deg2rad(10));
maxBlend = getdef_nested(params,'fb','max_blend',0.7);

% --- [SEATBELT TOGGLE] default ON, overridable by params and 4th arg
useSeatbelt = getdef_nested(params,'fb','enableSeatbelt', true);
if ~isempty(varargin)
    flag = varargin{1};
    if islogical(flag)
        useSeatbelt = flag;
    elseif ischar(flag) || isstring(flag)
        switch lower(string(flag))
            case {"on","seatbelt","enable","true"}
                useSeatbelt = true;
            case {"off","noseatbelt","disable","false"}
                useSeatbelt = false;
            otherwise
                warning('M_block: unknown seatbelt flag ''%s''; keeping default.', string(flag));
        end
    else
        warning('M_block: unsupported seatbelt flag type; keeping default.');
    end
end

% ---------- persistent ghost state, smoothers & stats
persistent gSt aN_s_prev jN_s_prev initd stats
if isempty(initd) || getdef(params,'resetFW',false)
    gSt = getdef(params,'fw_init_state',zeros(13,1));    % [x y z u v w q0 q1 q2 q3 p q r]'
    aN_s_prev=zeros(3,1); jN_s_prev=zeros(3,1);
    stats = struct('quatReset',0,'nanClamped',0,'illCondSolve',0,'tauSaturated',0,'gainShapeFix',0,'demandClamped',0);
    initd = true;
end

% ---------- advance ghost FW (RK4) with same u_fw, then normalize
u = u_fw(:).';
k1 = fw_6dof_quat(gSt,          u(1), u(2), u(3), u(4), params);
k2 = fw_6dof_quat(gSt+0.5*dt*k1,u(1), u(2), u(3), u(4), params);
k3 = fw_6dof_quat(gSt+0.5*dt*k2,u(1), u(2), u(3), u(4), params);
k4 = fw_6dof_quat(gSt+dt*k3,    u(1), u(2), u(3), u(4), params);
gSt = gSt + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
% normalize quaternion
qg = gSt(7:10); nq = norm(qg);
if ~isfinite(nq) || nq < 1e-9
    gSt(7:10) = [1;0;0;0]; stats.quatReset = stats.quatReset + 1;
else
    gSt(7:10) = qg/nq;
end

% derivatives at post-step for refs (a*, ω̇*)
f_g = fw_6dof_quat(gSt, u(1), u(2), u(3), u(4), params);

% ---------- ghost refs
vb_g       = gSt(4:6);               % body vel
qBN_g      = gSt(7:10);
omega_g    = gSt(11:13);
vb_dot_g   = f_g(4:6);
omegadot_g = f_g(11:13);

R_BN_g = quatNED_to_body(qBN_g);  R_NB_g = R_BN_g.';
aN_g_raw = R_NB_g * ( vb_dot_g + cross(omega_g, vb_g) );

% IIR smoothing on a* then j*, then s*
aN_s  = smooth_IIR(aN_g_raw, aN_s_prev, alpha_a);
jN_s  = smooth_IIR((aN_s - aN_s_prev)/dt, jN_s_prev, alpha_j);
sN_g  = (jN_s - jN_s_prev)/dt;
aN_s_prev = aN_s;  jN_s_prev = jN_s;

psi_g    = atan2(R_NB_g(2,1), R_NB_g(1,1));
psidot_g = [0 0 1] * (R_NB_g * omega_g);
psiddot_g= [0 0 1] * (R_NB_g * omegadot_g);

% ---------- QC current geometry & kinematics
R_BN = quatNED_to_body(qBN_qc);  R_NB = R_BN.';
e3 = [0;0;1];  b3 = R_NB(:,3);
omega = [p;q;r];

b3dot = R_NB * cross(omega, e3);
Jomega = J*omega;
omegadot_bias = J \ (-cross(omega, Jomega));

% b3ddot split = bias + (tau-part)
b3ddot_bias = R_NB * ( cross(omegadot_bias, e3) + cross(omega, cross(omega, e3)) );
M = [0 1 0; -1 0 0; 0 0 0];                       % a × e3 = M a
b3ddot_tau  = R_NB * ( M * (J \ eye(3)) );

% model-based a, j for QC (no drag)
aN_qc = -(T/m)*b3 + g*e3;
jN_qc = -(nu_T/m)*b3 - (T/m)*b3dot;

% ---------- outer-loop shaping (4th order)
e_r = [gSt(1);gSt(2);gSt(3)] - [xq;yq;zq];
e_v = (R_NB_g*vb_g) - [vN;vE;vD];
e_a = aN_s - aN_qc;
e_j = jN_s - jN_qc;
sN_cmd = sN_g + gmul(k3,e_j) + gmul(k2,e_a) + gmul(k1,e_v) + gmul(k0,e_r);
[sN_cmd, wasC1] = cap_vec(sN_cmd, s_max);  % snap cap

% --- thrust-feasible accel clamp on a* (keep ||g e3 - a*|| in [Tmin/m, Tmax/m])
Tmin     = getdef(params,'qc_T_min',0);
Tmax     = getdef(params,'qc_T_max',inf);
v_des    = g*e3 - aN_s;                             % proportional to thrust/m
nv       = norm(v_des);
nv_clamp = min(max(nv, Tmin/max(m,1e-9)), Tmax/max(m,1e-9));
if isfinite(nv_clamp) && nv > 0
    v_des = v_des * (nv_clamp / nv);
    aN_s  = g*e3 - v_des;                           % feasible a*
end

% --- authority/yaw measures based on desired force direction
f_des     = m*( g*e3 - aN_s );                      % desired force (NED, Z-down)
nf        = norm(f_des);
cos_gamma = (f_des.' * b3) / max(nf,1e-9);          % alignment b3 ↔ desired force
mg        = m*g;

% yaw 2nd-order shaping (ψ̈ command) with authority-aware cap
psi_qc    = atan2(R_NB(2,1), R_NB(1,1));
psidot_qc = [0 0 1] * (R_NB * omega);
e_psi  = wrapToPi(psi_g - psi_qc);
e_psid = psidot_g - psidot_qc;
psidd_cmd = psiddot_g + scalar(kpsi1)*e_psid + scalar(kpsi0)*e_psi;

% yaw authority weighting (0.3..1)
wpsi = max(0.3, min([ b3(3), cos_gamma, T/max(mg,1e-9) ]));
psidd_cmd = cap_scalar(psidd_cmd, wpsi * psidd_max);

% ---------- Safety attitude feedback (QC -> ghost), blended 0..maxBlend
R_des = R_NB_g;          % desired orientation (ghost FW)
R     = R_NB;
Re    = R_des.' * R;
e_R   = 0.5 * [ Re(3,2)-Re(2,3); Re(1,3)-Re(3,1); Re(2,1)-Re(1,2) ];  % SO(3) error
omega_d = omega_g;                                                     % ghost body rates
e_Om  = omega - (R.' * R_des) * omega_d;

theta_err = norm(e_R);  % keep for diagnostics
if useSeatbelt
    g1 = min(1, theta_err / max(thetaAct,1e-6));                       % from attitude error
    g2 = 1 - max(0.0, min([ b3(3), T/max(mg,1e-9), cos_gamma ]));      % from low authority
    gamma_fb = maxBlend * max(g1, g2);                                  % 0..maxBlend
else
    gamma_fb = 0;                                                       % seatbelt off
end

% [SEATBELT] Anti-windup throttling (only when enabled)
if useSeatbelt
    sN_cmd    = (1 - 0.5*gamma_fb) * sN_cmd;
    psidd_cmd = (1 - 0.5*gamma_fb) * psidd_cmd;
end

% ---------- 4×4 decoupling solve
alpha_top = -(2*nu_T/m)*b3dot - (T/m)*b3ddot_bias;                      % 3x1
alpha_yaw = [0 0 1] * (R_NB * (J \ (-cross(omega, Jomega))));           % scalar

B_top = [ -(1/m)*b3,  -(T/m)*b3ddot_tau ];                              % 3×(1+3)
aPsi  = [0 0 1] * R_NB * (J \ eye(3));                                  % 1×3

% Row normalization + yaw weighting → keeps τz from spiking when authority is tight
ntop = norm(B_top,'fro')/sqrt(3);                                       % avg top-row norm
nyaw = norm(aPsi);
row_scale = ntop / max(nyaw,1e-9);

B = [ B_top;
      0, (row_scale*wpsi) * aPsi ];

d = [ sN_cmd - alpha_top;
      (row_scale*wpsi) * (psidd_cmd - alpha_yaw) ];

% diag-scaled ridge
colN = sqrt(sum(B.^2,1))';  Dlam = diag(max(colN,1e-6));
b3z = b3(3);

% smooth |T| used in lambda_T term (prevents λ jumps)
persistent T_for_lam
if isempty(T_for_lam), T_for_lam = abs(T); end
alpha_lam = 0.2;
T_for_lam = (1-alpha_lam)*T_for_lam + alpha_lam*abs(T);

lambda = lam0 + lam_b3*(1 - b3z^2) + lam_T*( Tref/(Tref + max(1e-6, T_for_lam)) );

BtB = B.'*B + lambda*(Dlam.'*Dlam);
rhs = B.'*d;
if ~all(isfinite(BtB),'all') || rcond(BtB) < 1e-9
    stats.illCondSolve = stats.illCondSolve + 1;
    u = pinv(B, 1e-3) * d;
else
    u = BtB \ rhs;
end
u = finite_or_zero(u);

Tdd_dec = u(1);
tau     = u(2:4);

% ---------- thrust feed-forward around T*  (vertical balance + magnitude)
T_star_mag  = norm(f_des);                                   % magnitude target
b3z_eff     = max(b3(3), 0.15);                              % guard tiny b3z
a_z_des     = aN_s(3);                                       % desired vertical accel (down+)
T_star_vert = m * (g - a_z_des) / b3z_eff;                   % vertical support target

T_star_raw  = max(T_star_vert, T_star_mag);
T_floor_frac= getdef_nested(params,'dfl','T_floor_frac',0.10);
T_min       = getdef(params,'qc_T_min',0);
T_max       = getdef(params,'qc_T_max',inf);
T_floor     = T_floor_frac * m * g;
T_star      = min( max(T_star_raw, max(T_floor, T_min)), T_max );

Tdd_ff    = kT0*(T_star - T) - kT1*nu_T;
Tddot_cmd = Tdd_dec + Tdd_ff;
if isfinite(Tdd_lim), Tddot_cmd = min(max(Tddot_cmd, -Tdd_lim), Tdd_lim); end
Tddot_cmd = finite_or_zero(Tddot_cmd);

% ---------- add geometric PD safety torque (only if enabled), then saturate
if useSeatbelt
    tau_corr = -(kR(:).*e_R) - (kOm(:).*e_Om);
    tau      = tau + gamma_fb * tau_corr;
end
tau = min(max(tau, -abs(tau_lim)), abs(tau_lim));
if any(abs(tau) >= abs(tau_lim) - 1e-12)
    stats.tauSaturated = stats.tauSaturated + 1;
end
tau_roll  = tau(1);
tau_pitch = tau(2);
tau_yaw   = tau(3);

% ---------- Gimbal POV match (R_CB = R_x(η)R_y(θg))
R_CB_des = R_BN_g * R_NB;                                   % desired camera->body
theta_des = atan2( R_CB_des(1,3), R_CB_des(1,1) );          % θg*
eta_des   = atan2( R_CB_des(3,2), R_CB_des(2,2) );          % η*
gmax   = getdef(params,'gim_ang_max',deg2rad([60;45]));
margin = deg2rad(5);
theta_des = cap_scalar(theta_des, gmax(2)-margin);
eta_des   = cap_scalar(eta_des,   gmax(1)-margin);

gimRollAccCmd  = kg0_r * wrapToPi(eta_des   - eta)  + kg1_r * (0 - etad);
gimPitchAccCmd = kg0_p * wrapToPi(theta_des - thg)  + kg1_p * (0 - thgd);

% ---------- diagnostics
Mskew = [0 1 0; -1 0 0; 0 0 0];
BetaTau  = (T/m) * (R_NB * (Mskew * (J \ eye(3))));         % τ→snap sensitivity
normBetaTau = norm(BetaTau,2);
normBTcol1  = norm(-(1/m)*b3);                              % T̈→snap sensitivity
rcondBtB    = rcond(BtB);                                   % solve conditioning

ghost_diag = struct( ...
    'T',T, 'nu_T',nu_T, ...
    'T_star',T_star, 'T_star_raw',T_star_raw, ...
    'T_star_vert',T_star_vert, 'T_star_mag',T_star_mag, ...
    'T_floor',T_floor, ...
    'b3z',b3(3), 'b3z_eff',b3z_eff, ...
    'lambda',lambda, ...
    'normBetaTau',normBetaTau, 'normBT',normBTcol1, ...
    'rcondBtB',rcondBtB, ...
    'wpsi',wpsi, 'cos_gamma',cos_gamma, ...
    'theta_err',theta_err, 'gamma_fb',gamma_fb, ...
    'useSeatbelt',logical(useSeatbelt) );

stats.demandClamped = stats.demandClamped + double(wasC1 || abs(psidd_cmd)>=psidd_max-1e-12);

% ---------- ghost bundle
ghost = struct();
ghost.rN   = [gSt(1); gSt(2); gSt(3)];
ghost.vN   = R_NB_g*vb_g;
ghost.aN   = aN_s;
ghost.jN   = jN_s;
ghost.sN   = sN_g;
ghost.psi     = psi_g;
ghost.psidot  = psidot_g;
ghost.psiddot = psiddot_g;
ghost.R_NB_fw  = R_NB_g;
ghost.R_CB_des = R_CB_des;
ghost.flags = stats;
ghost.diag  = ghost_diag;

% ================== nested helpers ==================
    function y = gmul(K, v)
        % Gain multiply with shape coercion and flagging.
        v = v(:);
        if isscalar(K), y = K*v; return; end
        if isequal(size(K),[3,3]), y = K*v; return; end
        if isvector(K)
            kvec = K(:);
            if numel(kvec) ~= 3
                stats.gainShapeFix = stats.gainShapeFix + 1;
                if numel(kvec) >= 3, kvec = kvec(1:3);
                else, kvec = [kvec; repmat(kvec(end), 3-numel(kvec), 1)]; end
            end
            y = kvec .* v; return;
        end
        stats.gainShapeFix = stats.gainShapeFix + 1;
        y = K(1:3,1:3) * v;
    end

    function y = finite_or_zero(x)
        y = x;
        bad = ~isfinite(y);
        if any(bad,'all')
            y(bad) = 0;
            stats.nanClamped = stats.nanClamped + 1;
        end
    end

    function y = smooth_IIR(curr, prev, alpha)
        y = (1-alpha)*prev + alpha*curr;
        if any(~isfinite(y))
            y = prev; % hold last good value
            stats.nanClamped = stats.nanClamped + 1;
        end
    end

    function [v,wasC] = cap_vec(v, vmax)
        vmax = vmax(:); v = v(:);
        if numel(vmax) == 1, vmax = repmat(vmax,3,1); end
        wasC = any(abs(v) > vmax);
        v = sign(v) .* min(abs(v), vmax);
    end

    function s = cap_scalar(s, smax)
        s = min(max(s, -smax), smax);
    end

    function s = scalar(x)
        if isscalar(x), s = x; else, s = x(1); end
    end
end

% ========================= file-scope helpers =========================
function R = quatNED_to_body(q)
% q = [q0 q1 q2 q3], Hamilton, N->B
q0=q(1); q1=q(2); q2=q(3); q3=q(4);
R = [ 1-2*(q2^2+q3^2),   2*(q1*q2 + q0*q3),  2*(q1*q3 - q0*q2);
      2*(q1*q2 - q0*q3), 1-2*(q1^2+q3^2),    2*(q2*q3 + q0*q1);
      2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1),  1-2*(q1^2+q2^2) ];
end

function out = getdef(S, field, defaultVal)
if isfield(S, field), out = S.(field); else, out = defaultVal; end
end

function out = getdef_nested(S, f1, f2, defaultVal)
if isfield(S,f1) && isfield(S.(f1),f2), out = S.(f1).(f2); else, out = defaultVal; end
end

function ang = wrapToPi(ang)
ang = mod(ang + pi, 2*pi) - pi;
end
