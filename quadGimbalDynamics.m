function xdot = quadGimbalDynamics(state, Tddot_cmd, tau_roll, tau_pitch, tau_yaw, ...
                                   gimRollAccCmd, gimPitchAccCmd, params)
% QUADGIMBALDYNAMICS  6-DoF quadrotor with thrust double-integrator + 2-axis gimbal
%
% State (19x1):
%   [x y z vN vE vD q0 q1 q2 q3 p q r T nu_T eta theta_g eta_dot theta_g_dot]'
% Inputs:
%   Tddot_cmd        : commanded T̈  (N/s^2)
%   tau_roll/pitch/yaw: body torques (N·m) about Bx,By,Bz
%   gimRollAccCmd    : commanded gimbal roll  acceleration η̈ (rad/s^2)
%   gimPitchAccCmd   : commanded gimbal pitch acceleration θ̈g (rad/s^2)
% Params (required fields; sensible defaults used if missing):
%   qc_m, qc_J, qc_T_min, qc_T_max, qc_tau_max (3x1), g, rho, dt
%   wind_ned (3x1), qc_CdA (1x3, NED axes)  [optional; default zeros]
%   gim_ang_max (2x1), gim_rate_max (2x1), gimbal_ideal (bool) [optional]
%
% Frames: NED (Z down). Total thrust acts along -z_B.
% Quaternion convention: scalar-first Hamilton, q_BN rotates NED→Body.

% -------- Unpack state
x     = state(1);  %#ok<NASGU>
y     = state(2);  %#ok<NASGU>
z     = state(3);  %#ok<NASGU>
vN    = state(4);
vE    = state(5);
vD    = state(6);
qBN   = state(7:10);     % [q0 q1 q2 q3], N->B
pqr   = state(11:13);    % [p q r] (rad/s), body rates
T     = state(14);       % total thrust (N)
nu_T  = state(15);       % Tdot (N/s)
eta   = state(16);       % gimbal roll about Bx (rad)
theta = state(17);       % gimbal pitch about By after Rx (rad)
etad  = state(18);       % η̇
thetad= state(19);       % θ̇g

% -------- Params & defaults
m   = getfield_def(params,'qc_m',1.0);
J   = getfield_def(params,'qc_J',diag([0.02 0.02 0.035]));
g   = getfield_def(params,'g',9.81);
rho = getfield_def(params,'rho',1.225);

Tmin = getfield_def(params,'qc_T_min',0);
Tmax = getfield_def(params,'qc_T_max',30);
tau_lim = getfield_def(params,'qc_tau_max',[1.2;1.2;1.2]);

CdA = getfield_def(params,'qc_CdA',[0 0 0]);     % per-axis (NED)
wind_ned = getfield_def(params,'wind_ned',[0;0;0]);

omega_damp = getfield_def(params,'qc_omega_damp',[0.02 0.02 0.02]);  % N·m per rad/s
D = diag(omega_damp);

% (Optional) thrust rate/acc limits
Tdot_max  = getfield_def(params,'qc_Tdot_max',inf);
Tddot_max = getfield_def(params,'qc_Tddot_max',inf);

% Gimbal limits
gim_ang_max  = getfield_def(params,'gim_ang_max',deg2rad([60;45]));   % [roll; pitch]
gim_rate_max = getfield_def(params,'gim_rate_max',deg2rad([360;360]));
gim_acc_max  = getfield_def(params,'gim_acc_max',inf(2,1));           % optional
gimbal_ideal = getfield_def(params,'gimbal_ideal',false);

% -------- Rotation stuff
R_BN = quatNED_to_body(qBN);   % N -> B
R_NB = R_BN.';                 % B -> N
b3_N = R_NB(:,3);              % body z-axis expressed in NED

% -------- Kinematics (positions)
x_dot = vN;
y_dot = vE;
z_dot = vD;

% -------- Translational dynamics (NED)
vNED   = [vN; vE; vD];
v_air  = vNED - wind_ned;

% Thrust in NED: F_thrust = -T * b3_N
F_thrust_N = -T * b3_N;

% Quadratic drag per-axis in NED (optional; zero by default)
F_drag_N = -0.5 * rho * [CdA(1)*abs(v_air(1))*v_air(1);
                         CdA(2)*abs(v_air(2))*v_air(2);
                         CdA(3)*abs(v_air(3))*v_air(3)];

aNED = (F_thrust_N + F_drag_N)/m + [0;0;g];   % g is +down in NED

vN_dot = aNED(1);
vE_dot = aNED(2);
vD_dot = aNED(3);

% -------- Attitude kinematics (quaternion) & rotational dynamics
p = pqr(1); q = pqr(2); r = pqr(3);
q0 = qBN(1); q1 = qBN(2); q2 = qBN(3); q3 = qBN(4);

Omega = [ 0   -p   -q   -r;
          p    0    r   -q;
          q   -r    0    p;
          r    q   -p    0 ];
qBN_dot = 0.5 * Omega * [q0;q1;q2;q3];
% (Normalize q in the driver after integration.)

% Body torques with saturation
tau_cmd = [tau_roll; tau_pitch; tau_yaw];
tau     = saturate_vec(tau_cmd, tau_lim);

% Euler's equation
pqr_dot = J \ ( tau - cross(pqr, J*pqr) - D*pqr );

% -------- Thrust double integrator with projection anti-windup
Tdd = saturate_scalar(Tddot_cmd, -Tddot_max, Tddot_max);

% Default (inside bounds)
T_dot   = nu_T;
nu_Tdot = Tdd;

% Hard projection at bounds
if T <= Tmin
    if nu_T <= 0,   T_dot = 0;        end           % don't integrate below Tmin
    if Tdd  <  0,   nu_Tdot = 0;      end           % don't accelerate further down
elseif T >= Tmax
    if nu_T >= 0,   T_dot = 0;        end           % don't integrate above Tmax
    if Tdd  >  0,   nu_Tdot = 0;      end
end


% -------- Gimbal joint dynamics (roll=η about Bx, pitch=θ about By after Rx)
eta_dd   = saturate_scalar(gimRollAccCmd,  -gim_acc_max(1),  gim_acc_max(1));
theta_dd = saturate_scalar(gimPitchAccCmd, -gim_acc_max(2),  gim_acc_max(2));

if ~gimbal_ideal
    % Rate limits (if at rate limit and trying to speed up further, block)
    if abs(etad)   >= gim_rate_max(1) && sign(eta_dd)   == sign(etad),   eta_dd   = 0; end
    if abs(thetad) >= gim_rate_max(2) && sign(theta_dd) == sign(thetad), theta_dd = 0; end

    % Angle soft stops: if beyond angle limit and still pushing outwards, damp it
    k_stop = 50;   % damping toward boundary when saturated
    if abs(eta) >= gim_ang_max(1) && sign(etad) == sign(eta)
        if sign(eta_dd) == sign(eta), eta_dd = -k_stop*etad; end
    end
    if abs(theta) >= gim_ang_max(2) && sign(thetad) == sign(theta)
        if sign(theta_dd) == sign(theta), theta_dd = -k_stop*thetad; end
    end
else
    % "Ideal" path: no extra limiting here (use very large limits in params)
end

% -------- Assemble xdot
xdot = zeros(19,1);
xdot(1)  = x_dot;
xdot(2)  = y_dot;
xdot(3)  = z_dot;
xdot(4)  = vN_dot;
xdot(5)  = vE_dot;
xdot(6)  = vD_dot;
xdot(7:10)  = qBN_dot;
xdot(11:13) = pqr_dot;
xdot(14) = T_dot;
xdot(15) = nu_Tdot;
xdot(16) = etad;
xdot(17) = thetad;
xdot(18) = eta_dd;
xdot(19) = theta_dd;

end

% ===================== helpers ==========================
function R = quatNED_to_body(q)
% q = [q0 q1 q2 q3] (scalar-first), mapping NED -> Body
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
% Hamilton, passive, N->B
R = [ 1-2*(q2^2+q3^2),   2*(q1*q2 + q0*q3),  2*(q1*q3 - q0*q2);
      2*(q1*q2 - q0*q3), 1-2*(q1^2+q3^2),    2*(q2*q3 + q0*q1);
      2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1),  1-2*(q1^2+q2^2) ];
end

function y = saturate_scalar(u, umin, umax)
y = min(max(u, umin), umax);
end

function y = saturate_vec(u, umax_vec)
% symmetric limits per axis: [-umax, +umax]
y = min(max(u, -abs(umax_vec)), abs(umax_vec));
end

function out = getfield_def(S, field, defaultVal)
if isfield(S, field), out = S.(field); else, out = defaultVal; end
end
