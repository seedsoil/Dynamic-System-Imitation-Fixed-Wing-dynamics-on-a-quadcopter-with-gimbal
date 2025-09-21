function S = run_fw_qc_compare(seq, opts)
% RUN_FW_QC_COMPARE  FW vs QC compare; pass 'seq' (Nx2 [channel, signedMag]).
%   S = run_fw_qc_compare(seq)
%   S = run_fw_qc_compare(seq, Name=Value, ...)
%
% Name-Value options
%   PlotsRoot        : base folder for exports (default: fullfile(pwd,'Plots'))
%   ExportVariants   : ["both","on"] or ["on"] (duplicates policy)
%   SavePNG, SavePDF : logical (true, true)
%   Verbose          : logical (true)   % also saved into evaluation.txt
%   MarkImpulses     : logical (false)  % default OFF: no vertical dotted lines
%   TsimTotal        : total sim time [s]           (30)
%   TAnalyzeTo       : analysis window end [s]      (20)
%   PulseWidth       : pulse width Tw [s]           (1.00)
%   FirstPulseCenter : first pulse center [s]       (2)
%   AutoFitSchedule  : spread pulses to end window  (true)
%   SettleMargin     : quiet time before end [s]    (2)
%   ThrustBias       : baseline thrust [N]          (0.4)
%   SurfBias         : baseline surface [rad]       (0)
%   AmpSurfacePerUnitDeg : fixed peak per unit [deg] for surfaces (2.0)
%   AmpThrustPerUnitN    : fixed peak per unit [N]   for thrust   (6.0)
%   ShowLimitLines   : logical (false)  % hides T-min/T-max lines unless asked
%   LegendsOutside   : logical (true)   % keep legends off the data

arguments
    seq (:,2) double
    opts.TsimTotal        (1,1) double  = 30.0
    opts.TAnalyzeTo       (1,1) double  = 20.0
    opts.PulseWidth       (1,1) double  = 1.00
    opts.FirstPulseCenter (1,1) double  = 2.00
    opts.AutoFitSchedule  (1,1) logical = true
    opts.SettleMargin     (1,1) double  = 2.0
    opts.ThrustBias       (1,1) double  = 0.4
    opts.SurfBias         (1,1) double  = 0.0
    opts.AmpSurfacePerUnitDeg (1,1) double = 2.0
    opts.AmpThrustPerUnitN    (1,1) double = 6.0
    opts.ExportVariants   (1,:) string  {mustBeMember(opts.ExportVariants,["both","on"])} = ["both","on"]
    opts.PlotsRoot        (1,1) string  = ""
    opts.SavePNG          (1,1) logical = true
    opts.SavePDF          (1,1) logical = true
    opts.MarkImpulses     (1,1) logical = false
    opts.Verbose          (1,1) logical = true
    opts.ShowLimitLines   (1,1) logical = false
    opts.LegendsOutside   (1,1) logical = true
end

% ---------- params from script (isolated wrapper) ----------
[params, initial_fw_state, initial_qc_state] = local_load_params();  % runs param.m safely
set(groot,'defaultLegendAutoUpdate','off');

% ---------- colors & style ----------
C_FW  = [0.0000 0.4470 0.7410];   % blue
C_OFF = [0.0000 0.7000 0.2000];   % greener, more vibrant
C_ON  = [0.8500 0.3250 0.0980];   % red
LW    = 1.3;

% ---------- global plot look: subtle, background grid ----------
set(groot, ...
    'defaultAxesXGrid','on', 'defaultAxesYGrid','on', ...                  % major grid on
    'defaultAxesGridLineStyle','-', ...                                     % solid grid
    'defaultAxesGridColor',[0.85 0.85 0.85], ...                            % light gray
    'defaultAxesGridAlpha',0.10, ...                                        % faint
    'defaultAxesXMinorGrid','off','defaultAxesYMinorGrid','off', ...        % no minor grid
    'defaultAxesLayer','bottom');                                           % grid behind data  % Axes Layer property

% ---------- time base ----------
dt = params.dt;
N  = floor(opts.TsimTotal/dt) + 1;
t  = (0:N-1)'*dt;
kinfo = t <= opts.TAnalyzeTo;
tP    = t(kinfo); nP = numel(tP);

% ---------- pulse width (Tw) and fixed peak mapping ----------
Tw = opts.PulseWidth;
A_surface_unit_rad = deg2rad(opts.AmpSurfacePerUnitDeg);
A_thrust_unit_N    = opts.AmpThrustPerUnitN;

% ---------- impulse schedule (centers) ----------
numP = size(seq,1);
if opts.AutoFitSchedule
    t0s = linspace(opts.FirstPulseCenter, opts.TAnalyzeTo - opts.SettleMargin, numP);
else
    t0s = opts.FirstPulseCenter + (0:numP-1)*2.0;
end
rect = @(ti,t0) double( (ti >= (t0 - 0.5*Tw)) & (ti < (t0 + 0.5*Tw)) );

% ---------- pulse heights: fixed per-unit peaks (independent of Tw) ----------
pulseHeights = zeros(numP,1);
for i=1:numP
    ch  = seq(i,1);
    k   = seq(i,2);
    if ch==1
        pulseHeights(i) = k * A_thrust_unit_N;          % [N]
    else
        pulseHeights(i) = k * A_surface_unit_rad;       % [rad]
    end
end

% ---------- ulog (FW command timeline) ----------
THRUST_BIAS_N = opts.ThrustBias;
SURF_BIAS_RAD = opts.SurfBias;
ulog = zeros(N,4);                             % [thrust elev ailer rud]
for k=1:N
    ti = t(k);
    u = [THRUST_BIAS_N; SURF_BIAS_RAD; 0; 0];
    if ti <= opts.TAnalyzeTo
        for i=1:numP
            amp = pulseHeights(i) * rect(ti, t0s(i));
            if amp~=0
                switch seq(i,1)
                    case 1, u(1)=u(1)+amp;
                    case 2, u(2)=u(2)+amp;
                    case 3, u(3)=u(3)+amp;
                    case 4, u(4)=u(4)+amp;
                end
            end
        end
    end
    ulog(k,:) = u(:).';
end

% ---------- FW trim & run ----------
Va_star = 15; alpha_trim = -0.034579507; elev_trim = 0.072040640; T_trim = 0.256063212;
THRUST_BIAS_N = T_trim; %#ok<NASGU>
ELEV_BIAS_RAD = elev_trim; %#ok<NASGU>
u0 = Va_star*cos(alpha_trim); w0 = Va_star*sin(alpha_trim);
initial_fw_state(4) = u0; initial_fw_state(5) = 0; initial_fw_state(6) = w0;
phi0=0; theta0=alpha_trim; psi0=0;
cph=cos(phi0/2); sph=sin(phi0/2); cth=cos(theta0/2); sth=sin(theta0/2); cps=cos(psi0/2); sps=sin(psi0/2);
q0=cph*cth*cps + sph*sth*sps; q1=cph*cth*sps - sph*sth*cps; q2=cph*sth*cps + sph*cth*sps; q3=sph*cth*cps - cph*sth*sps;
initial_fw_state(7:10)=[q0;q1;q2;q3]/norm([q0;q1;q2;q3]); initial_fw_state(11:13)=[0;0;0];

Xfw = zeros(N,13); Xfw(1,:) = initial_fw_state(:).';
Va_fw=zeros(N,1); psi_fw=zeros(N,1); alt_fw=zeros(N,1);
phi_fw=zeros(N,1); theta_fw=zeros(N,1);
for k=1:N-1
    xf = Xfw(k,:).';  u = ulog(k,:).';
    k1 = fw_6dof_quat(xf,               u(1),u(2),u(3),u(4),params);
    k2 = fw_6dof_quat(xf + 0.5*dt*k1,   u(1),u(2),u(3),u(4),params);
    k3 = fw_6dof_quat(xf + 0.5*dt*k2,   u(1),u(2),u(3),u(4),params);
    k4 = fw_6dof_quat(xf + dt*k3,       u(1),u(2),u(3),u(4),params);
    xfn = xf + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    qf = xfn(7:10); nq=norm(qf); qf=qf/max(nq,1e-12); xfn(7:10)=qf;
    Xfw(k+1,:) = xfn.';
end
for k=1:N
    qBN_f = Xfw(k,7:10).'; R_BN_f = quatNED_to_body(qBN_f); R_NB_f = R_BN_f.';
    [phi_fw(k), theta_fw(k), psi_fw(k)] = euler321_from_RNB(R_NB_f);
    alt_fw(k) = Xfw(k,3);
    v_b  = Xfw(k,4:6).'; Va_fw(k) = norm(v_b - R_BN_f*params.wind_ned);
end
psi_fw = unwrap(psi_fw);

% ---------- run QC (off / on) ----------
flags = struct('watchdog',true,'stop_on_bad',true,'mark_impulses',opts.MarkImpulses, ...
               'debug_print',false,'authority_trip',true,'rcond_thresh',5e-6,'b3z_thresh',0.15, ...
               'T_floor_guard_frac',0.02,'on_trip_action','warn','clipTau_margin',0.05, ...
               'scaleTau_factor',0.50,'trip_once_only',false);
if ~isfield(params,'dfl'), params.dfl = struct(); end
if ~isfield(params.dfl,'T_floor_frac'), params.dfl.T_floor_frac = 0.10; end
if ~isfield(params.dfl,'lambda_T'),     params.dfl.lambda_T     = 0.6;  end
if ~isfield(params,'caps') || ~isfield(params.caps,'psidd_max'), params.caps.psidd_max = 0.50; end
if ~isfield(params,'caps') || ~isfield(params.caps,'s_max'),     params.caps.s_max     = [8;8;8];   end

qc_off = run_qc_once('off', params, initial_qc_state, ulog, t, dt, flags);
qc_on  = run_qc_once('on',  params, initial_qc_state, ulog, t, dt, flags);

% ---------- derived views ----------
[phi_off,theta_off,psi_off,phi_cam_off,theta_cam_off,psi_cam_off,pov_err_off] = derive_qc_views(qc_off.X, Xfw, params);
[phi_on ,theta_on ,psi_on ,phi_cam_on ,theta_cam_on ,psi_cam_on ,pov_err_on ] = derive_qc_views(qc_on.X , Xfw, params);
psi_off = unwrap(psi_off); psi_on = unwrap(psi_on);

% ---------- analysis windows ----------
XfwP = Xfw(kinfo,:);  XoffP = qc_off.X(kinfo,:);  XonP  = qc_on.X(kinfo,:);
phi_fwP=phi_fw(kinfo);  th_fwP=theta_fw(kinfo);  ps_fwP=psi_fw(kinfo);
phi_offP=phi_off(kinfo); th_offP=theta_off(kinfo); ps_offP=psi_off(kinfo);
phi_onP =phi_on(kinfo);  th_onP =theta_on(kinfo);  ps_onP =psi_on(kinfo);
pov_err_offP = pov_err_off(kinfo);  pov_err_onP = pov_err_on(kinfo);

% Altitude up-positive
alt_fw_up = -XfwP(:,3); alt_off_up = -XoffP(:,3); alt_on_up = -XonP(:,3);

% ---------- SAFE getters (robust) ----------
T_offP     = get_or(qc_off,'T',        kinfo, XoffP(:,14));
T_onP      = get_or(qc_on ,'T',        kinfo, XonP(:,14));
nuT_offP   = get_or(qc_off,'nuT',      kinfo, XoffP(:,15));
nuT_onP    = get_or(qc_on ,'nuT',      kinfo, XonP(:,15));

Tstar_offP = get_or(qc_off,'T_star',   kinfo, nan(size(tP)));
Tstar_onP  = get_or(qc_on ,'T_star',   kinfo, nan(size(tP)));
Tddot_offP = get_or(qc_off,'Tddot',    kinfo, nan(size(tP)));
Tddot_onP  = get_or(qc_on ,'Tddot',    kinfo, nan(size(tP)));

tau_offP   = get_or(qc_off,'tau',      kinfo, nan(nP,3));
tau_onP    = get_or(qc_on ,'tau',      kinfo, nan(nP,3));

b3z_offP   = get_or(qc_off,'b3z',      kinfo, nan(size(tP)));
b3z_onP    = get_or(qc_on ,'b3z',      kinfo, nan(size(tP)));
rcond_offP = get_or(qc_off,'rcondBtB', kinfo, nan(size(tP)));
rcond_onP  = get_or(qc_on ,'rcondBtB', kinfo, nan(size(tP)));
lam_offP   = get_or(qc_off,'lambda',   kinfo, nan(size(tP)));
lam_onP    = get_or(qc_on ,'lambda',   kinfo, nan(size(tP)));

gamma_onP    = get_or(qc_on ,'gamma',    kinfo, nan(size(tP)));
wpsi_onP     = get_or(qc_on ,'wpsi',     kinfo, nan(size(tP)));
g1_onP       = get_or(qc_on ,'g1',       kinfo, nan(size(tP)));
g2_onP       = get_or(qc_on ,'g2',       kinfo, nan(size(tP)));
thetaErr_onP = get_or(qc_on,'theta_err', kinfo, nan(size(tP)));
cosgam_onP   = get_or(qc_on,'cos_gamma', kinfo, nan(size(tP)));

% ---------- metrics ----------
dx_off = XoffP(:,1)-XfwP(:,1); dy_off = XoffP(:,2)-XfwP(:,2); dz_off = XoffP(:,3)-XfwP(:,3);
dx_on  = XonP(:,1) -XfwP(:,1); dy_on  = XonP(:,2) -XfwP(:,2); dz_on  = XonP(:,3) -XfwP(:,3);
dpos3_off = sqrt(dx_off.^2 + dy_off.^2 + dz_off.^2);
dpos3_on  = sqrt(dx_on .^2 + dy_on .^2 + dz_on .^2);
wrap = @(a) mod(a + pi, 2*pi) - pi;
MAE_att_off = rad2deg([mean(abs(wrap(phi_offP-phi_fwP))) mean(abs(wrap(th_offP-th_fwP))) mean(abs(wrap(ps_offP-ps_fwP)))]);
MAX_att_off = rad2deg([max(abs(wrap(phi_offP-phi_fwP)))  max(abs(wrap(th_offP-th_fwP)))  max(abs(wrap(ps_offP-ps_fwP))) ]);
MAE_att_on  = rad2deg([mean(abs(wrap(phi_onP-phi_fwP ))) mean(abs(wrap(th_onP-th_fwP ))) mean(abs(wrap(ps_onP-ps_fwP )))]);
MAX_att_on  = rad2deg([max(abs(wrap(phi_onP-phi_fwP )))  max(abs(wrap(th_onP-th_fwP )))  max(abs(wrap(ps_onP-ps_fwP ))) ]);
MAE3_off = mean(dpos3_off); MAX3_off = max(dpos3_off);
MAE3_on  = mean(dpos3_on ); MAX3_on  = max(dpos3_on );
MAE_pov_off = mean(pov_err_offP,'omitnan'); MAX_pov_off = max(pov_err_offP,[],'omitnan');
MAE_pov_on  = mean(pov_err_onP ,'omitnan'); MAX_pov_on  = max(pov_err_onP ,[],'omitnan');

imp = @(on,off) 100*(off-on)/max(off,eps);
imp_pos  = imp(MAE3_on , MAE3_off);
imp_alt  = imp(mean(abs(alt_on_up-alt_fw_up)), mean(abs(alt_off_up-alt_fw_up)));
imp_attR = imp(MAE_att_on(1), MAE_att_off(1));
imp_attP = imp(MAE_att_on(2), MAE_att_off(2));
imp_attY = imp(MAE_att_on(3), MAE_att_off(3));
imp_pov  = imp(MAE_pov_on , MAE_pov_off);

% ---------- output folder ----------
seqToken = seq_to_token(seq);
plotsRoot = char( ifelse(strlength(opts.PlotsRoot)==0, fullfile(pwd,'Plots'), opts.PlotsRoot) );
outDir = fullfile(plotsRoot, seqToken); if ~isfolder(outDir), mkdir(outDir); end

% ---------- evaluation text ----------
evalText = compose_eval_text(opts.TAnalyzeTo, MAE3_off,MAE_xyz(XoffP,XfwP),MAX3_off, ...
                                      MAE_att_off,MAX_att_off,MAE_pov_off,MAX_pov_off, ...
                                      MAE3_on ,MAE_xyz(XonP ,XfwP),MAX3_on , ...
                                      MAE_att_on ,MAX_att_on ,MAE_pov_on ,MAX_pov_on , ...
                                      imp_pos,imp_alt,imp_attR,imp_attP,imp_attY,imp_pov, seqToken);
if opts.Verbose, fprintf('%s', evalText); end
write_text(fullfile(outDir,'evaluation.txt'), evalText);

% ---------- defaults for figures/legends/fonts ----------
set(groot,'defaultAxesFontName','Helvetica','defaultAxesFontSize',10, ...
          'defaultLineLineWidth',LW,'defaultFigureColor','w');
set(groot,'defaultLegendFontSizeMode','manual','defaultLegendFontSize',10);

% ========== NEW: u_FW timeline (full figure; thrust vs surfaces) ==========
ufwP = ulog(kinfo,:);
f = figure('Name','u_{FW} timeline (commands)'); hold on;
yyaxis left
plot(tP, ufwP(:,1), 'Color',[0 0 0], 'DisplayName','T [N]'); ylabel('T [N]');
yyaxis right
plot(tP, rad2deg(ufwP(:,2)), '-', 'Color',[0.2 0.5 0.8], 'DisplayName','elev [deg]');
plot(tP, rad2deg(ufwP(:,3)), '-', 'Color',[0.8 0.4 0.4], 'DisplayName','ail [deg]');
plot(tP, rad2deg(ufwP(:,4)), '-', 'Color',[0.4 0.8 0.4], 'DisplayName','rud [deg]');
ylabel('[deg]'); grid on; xlabel('t [s]');
title('FW inputs u_{FW}: thrust (N) & surfaces (deg)'); legend('Location','best');
if opts.MarkImpulses, mark_impulses(t0s); end
saveplot(gcf, outDir, "u_FW_timeline", opts); close(f);

% POV axes — separate full figures (variants both/on)
plotPOVaxis(tP, rad2deg(phi_cam_off(kinfo)), rad2deg(phi_cam_on(kinfo)), rad2deg(phi_fwP), ...
    "POV Roll \phi_{cam} [deg]", "POV Roll", outDir, C_FW,C_OFF,C_ON, opts.ExportVariants, t0s, opts);
plotPOVaxis(tP, rad2deg(theta_cam_off(kinfo)), rad2deg(theta_cam_on(kinfo)), rad2deg(th_fwP), ...
    "POV Pitch \theta_{cam} [deg]", "POV Pitch", outDir, C_FW,C_OFF,C_ON, opts.ExportVariants, t0s, opts);
plotPOVaxis(tP, rad2deg(psi_cam_off(kinfo)), rad2deg(psi_cam_on(kinfo)), rad2deg(ps_fwP), ...
    "POV Yaw \psi_{cam} [deg]", "POV Yaw", outDir, C_FW,C_OFF,C_ON, opts.ExportVariants, t0s, opts);

% POV geodesic error — full figure
for variant = opts.ExportVariants
    f = figure('Name',"POV geodesic error — " + variant); hold on;
    if variant=="both"
        plot(tP, pov_err_offP,'--','Color',C_OFF);
        plot(tP, pov_err_onP , '-.','Color',C_ON); legend('QC off','QC on','Location','best');
    else
        plot(tP, pov_err_onP , '-','Color',C_ON);   legend('QC on','Location','best');
    end
    grid on; xlabel('t [s]'); ylabel('POV error [deg]'); title('POV geodesic error');
    if opts.MarkImpulses, mark_impulses(t0s); end
    saveplot(gcf, outDir, "POV_geodesic_error__" + variant, opts); close(f);
end

% NED position axes — separate full figures
plotNEDaxis(tP, XoffP(:,1), XonP(:,1), XfwP(:,1), 'x_N [m]', 'Position — North x_N', outDir, C_FW,C_OFF,C_ON, opts.ExportVariants, t0s, opts);
plotNEDaxis(tP, XoffP(:,2), XonP(:,2), XfwP(:,2), 'y_E [m]', 'Position — East y_E',  outDir, C_FW,C_OFF,C_ON, opts.ExportVariants, t0s, opts);
plotNEDaxis(tP, XoffP(:,3), XonP(:,3), XfwP(:,3), 'z_D [m] (down +)', 'Position — Down z_D', outDir, C_FW,C_OFF,C_ON, opts.ExportVariants, t0s, opts);

% Altitude (up+)
for variant = opts.ExportVariants
    f = figure('Name',"Altitude — " + variant); hold on;
    if variant=="both"
        plot(tP, alt_off_up,'--','Color',C_OFF);
        plot(tP, alt_on_up , '-.','Color',C_ON);
        plot(tP, alt_fw_up , '-','Color',C_FW);
        legend('QC off','QC on','FW','Location','best');
    else
        plot(tP, alt_on_up , '-','Color',C_ON);
        plot(tP, alt_fw_up , '-','Color',C_FW);
        legend('QC on','FW','Location','best');
    end
    grid on; xlabel('t [s]'); ylabel('Altitude [m] (up +)');
    title('Altitude (up-positive)'); if opts.MarkImpulses, mark_impulses(t0s); end
    saveplot(gcf, outDir, "Altitude_up__" + variant, opts); close(f);
end

% QC Thrust vs Target — full plot (both/on)
for variant = opts.ExportVariants
    f = figure('Name',"QC Thrust vs Target — " + variant); hold on;
    if variant=="both"
        plot(tP, T_offP,'--','Color',C_OFF);
        plot(tP, T_onP , '-.','Color',C_ON);
        if any(isfinite(Tstar_offP)), plot(tP, Tstar_offP,':','Color',C_OFF,'LineWidth',1.0); end
        if any(isfinite(Tstar_onP )), plot(tP, Tstar_onP ,':','Color',C_ON ,'LineWidth',1.0); end
        leg = {'T off','T on'}; if any(isfinite(Tstar_offP)), leg{end+1}='T* off'; end
        if any(isfinite(Tstar_onP)), leg{end+1}='T* on'; end
    else
        plot(tP, T_onP , '-','Color',C_ON);
        if any(isfinite(Tstar_onP)), plot(tP, Tstar_onP ,':','Color',C_ON,'LineWidth',1.0); end
        leg = {'T on'}; if any(isfinite(Tstar_onP)), leg{end+1}='T* on'; end
    end
    grid on; xlabel('t [s]'); ylabel('T [N]'); legend(leg,'Location','best');
    title('QC Thrust vs Target');

    % (Hide horizontal T_min/T_max lines unless requested)
    if opts.ShowLimitLines
        if isfield(params,'qc_T_min'), yline(params.qc_T_min, ':','HandleVisibility','off'); end
        if isfield(params,'qc_T_max'), yline(params.qc_T_max, ':','HandleVisibility','off'); end
    end

    if opts.MarkImpulses, mark_impulses(t0s); end
    saveplot(gcf, outDir, "QC_Thrust_vs_Target__" + variant, opts); close(f);
end

% AAGDA gamma(t) + Net added torque (ON-only)
f = figure('Name','AAGDA gamma(t) blend + Net added torque'); hold on;
yyaxis left;  plot(tP, gamma_onP, 'Color',C_ON); ylabel('\gamma(t) [–]');
yyaxis right; dtau = vecnorm(tau_onP - tau_offP,2,2);
plot(tP, dtau, 'Color',[0.25 0.25 0.25]); ylabel('|Δ\tau| [N m]');
xlabel('t [s]'); grid on; title('AAGDA \gamma(t) blend & Net Added Torque (ON run)');
if opts.MarkImpulses, mark_impulses(t0s); end
legend({'AAGDA \gamma(t)','|{\tau}_{on} - {\tau}_{off}|'},'Location','best');
saveplot(gcf, outDir, "AAGDA_gamma_and_NetAddedTorque__on", opts); close(f);

% Control parameters — Attitude (roll/pitch) — ON-only
f = figure('Name','Control Params — Attitude (roll/pitch)'); hold on;
yyaxis left
plot(tP, gamma_onP, 'Color',C_ON, 'DisplayName','\gamma(t)');
plot(tP, g1_onP,   '--','Color',[0.3 0.3 0.8],'DisplayName','g_1');
plot(tP, g2_onP,   ':', 'Color',[0.9 0.6 0.2],'DisplayName','g_2');
ylabel('[–]'); ylim([0 1.05]);
yyaxis right
plot(tP, rad2deg(thetaErr_onP), '-.','Color',[0.2 0.2 0.2],'DisplayName','\theta_{err} [deg]');
ylabel('[deg]');
grid on; xlabel('t [s]'); title('Control Params — Attitude (roll/pitch) — ON');
legend('Location','best'); if opts.MarkImpulses, mark_impulses(t0s); end
saveplot(gcf, outDir, "ControlParams_Attitude__on", opts); close(f);

% Control parameters — Yaw — ON-only
f = figure('Name','Control Params — Yaw'); hold on;
plot(tP, wpsi_onP,   '-', 'Color',C_ON,             'DisplayName','w_\psi');
plot(tP, cosgam_onP, '--','Color',[0.1 0.6 0.6],    'DisplayName','cos \gamma');
plot(tP, b3z_onP,    ':', 'Color',[0.6 0.4 0.1],    'DisplayName','b_{3z}');
grid on; xlabel('t [s]'); ylabel('[–]'); title('Control Params — Yaw — ON');
legend('Location','best'); if opts.MarkImpulses, mark_impulses(t0s); end
saveplot(gcf, outDir, "ControlParams_Yaw__on", opts); close(f);

%% === QC ON vs FW — POV & NED, 3x2 compact (no legends) ===
f = figure('Name','QC ON vs FW — POV & NED');
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% Row 1 — POV roll & POV pitch
nexttile; hold on;
plot(tP, rad2deg(phi_cam_on(kinfo)),'-.','Color',C_ON );
plot(tP, rad2deg(phi_fwP),          '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\phi_{cam} [deg]'); title('POV Roll');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, rad2deg(theta_cam_on(kinfo)),'-.','Color',C_ON );
plot(tP, rad2deg(th_fwP),             '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\theta_{cam} [deg]'); title('POV Pitch');
if opts.MarkImpulses, mark_impulses(t0s); end

% Row 2 — POV yaw & NED: North x_N
nexttile; hold on;
plot(tP, rad2deg(psi_cam_on(kinfo)),'-.','Color',C_ON );
plot(tP, rad2deg(ps_fwP),            '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\psi_{cam} [deg]'); title('POV Yaw');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, XonP(:,1),'-.','Color',C_ON );
plot(tP, XfwP(:,1),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('x_N [m]'); title('North x_N');
if opts.MarkImpulses, mark_impulses(t0s); end

% Row 3 — NED: East y_E & Down z_D
nexttile; hold on;
plot(tP, XonP(:,2),'-.','Color',C_ON );
plot(tP, XfwP(:,2),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('y_E [m]'); title('East y_E');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, XonP(:,3),'-.','Color',C_ON );
plot(tP, XfwP(:,3),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('z_D [m] (down +)'); title('Down z_D');
if opts.MarkImpulses, mark_impulses(t0s); end

title(tlo, 'QC ON vs FW — POV (roll,pitch,yaw) & NED (x,y,z) — 0..' + string(opts.TAnalyzeTo) + ' s');
saveplot(gcf, outDir, "QC_on_vs_FW__POV_NED__3x2_nolegends", opts);
close(f);

%% === QC OFF vs FW — POV & NED, 3x2 compact (no legends) ===
f = figure('Name','QC OFF vs FW — POV & NED');
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% Row 1 — POV roll & POV pitch
nexttile; hold on;
plot(tP, rad2deg(phi_cam_off(kinfo)),'--','Color',C_OFF);
plot(tP, rad2deg(phi_fwP),           '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\phi_{cam} [deg]'); title('POV Roll');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, rad2deg(theta_cam_off(kinfo)),'--','Color',C_OFF);
plot(tP, rad2deg(th_fwP),              '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\theta_{cam} [deg]'); title('POV Pitch');
if opts.MarkImpulses, mark_impulses(t0s); end

% Row 2 — POV yaw & NED: North x_N
nexttile; hold on;
plot(tP, rad2deg(psi_cam_off(kinfo)),'--','Color',C_OFF);
plot(tP, rad2deg(ps_fwP),             '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\psi_{cam} [deg]'); title('POV Yaw');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, XoffP(:,1),'--','Color',C_OFF);
plot(tP, XfwP(:,1),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('x_N [m]'); title('North x_N');
if opts.MarkImpulses, mark_impulses(t0s); end

% Row 3 — NED: East y_E & Down z_D
nexttile; hold on;
plot(tP, XoffP(:,2),'--','Color',C_OFF);
plot(tP, XfwP(:,2),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('y_E [m]'); title('East y_E');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, XoffP(:,3),'--','Color',C_OFF);
plot(tP, XfwP(:,3),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('z_D [m] (down +)'); title('Down z_D');
if opts.MarkImpulses, mark_impulses(t0s); end

title(tlo, 'QC OFF vs FW — POV (roll,pitch,yaw) & NED (x,y,z) — 0..' + string(opts.TAnalyzeTo) + ' s');
saveplot(gcf, outDir, "QC_off_vs_FW__POV_NED__3x2_nolegends", opts);
close(f);

%% === QC vs FW — POV & NED (off+on), 3x2 compact (no legends) ===
f = figure('Name','QC vs FW — POV & NED (off/on)');
tlo = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% Row 1 — POV roll & POV pitch
nexttile; hold on;
plot(tP, rad2deg(psi_cam_off(kinfo)),'--','Color',C_OFF);   % (typo fix not needed, left as in your script)
plot(tP, rad2deg(phi_cam_on (kinfo)),'-.','Color',C_ON );
plot(tP, rad2deg(phi_fwP),          '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\phi_{cam} [deg]'); title('POV Roll');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, rad2deg(theta_cam_off(kinfo)),'--','Color',C_OFF);
plot(tP, rad2deg(theta_cam_on (kinfo)),'-.','Color',C_ON );
plot(tP, rad2deg(th_fwP),             '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\theta_{cam} [deg]'); title('POV Pitch');
if opts.MarkImpulses, mark_impulses(t0s); end

% Row 2 — POV yaw & NED: North x_N
nexttile; hold on;
plot(tP, rad2deg(psi_cam_off(kinfo)),'--','Color',C_OFF);
plot(tP, rad2deg(psi_cam_on (kinfo)),'-.','Color',C_ON );
plot(tP, rad2deg(ps_fwP),            '-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('\psi_{cam} [deg]'); title('POV Yaw');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, XoffP(:,1),'--','Color',C_OFF);
plot(tP, XonP (:,1),'-.','Color',C_ON );
plot(tP, XfwP (:,1),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('x_N [m]'); title('North x_N');
if opts.MarkImpulses, mark_impulses(t0s); end

% Row 3 — NED: East y_E & Down z_D
nexttile; hold on;
plot(tP, XoffP(:,2),'--','Color',C_OFF);
plot(tP, XonP (:,2),'-.','Color',C_ON );
plot(tP, XfwP (:,2),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('y_E [m]'); title('East y_E');
if opts.MarkImpulses, mark_impulses(t0s); end

nexttile; hold on;
plot(tP, XoffP(:,3),'--','Color',C_OFF);
plot(tP, XonP (:,3),'-.','Color',C_ON );
plot(tP, XfwP (:,3),'-', 'Color',C_FW );
grid on; xlabel('t [s]'); ylabel('z_D [m] (down +)'); title('Down z_D');
if opts.MarkImpulses, mark_impulses(t0s); end

title(tlo, 'QC vs FW — POV (roll,pitch,yaw) & NED (x,y,z) — 0..' + string(opts.TAnalyzeTo) + ' s');
saveplot(gcf, outDir, "QC_vs_FW__POV_NED__off_on__3x2_nolegends", opts);
close(f);

% ========== NEW: u_QC timeline (M_block → QC dynamics), ON-only ==========
etaacc_onP = get_or(qc_on,'etaacc',kinfo, nan(size(tP)));
theacc_onP = get_or(qc_on,'theacc',kinfo, nan(size(tP)));
f = figure('Name','u_{QC} timeline (M\_block \rightarrow dynamics) — ON');
tlo2 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% top tile
nexttile; hold on;
yyaxis left
plot(tP, tau_onP(:,1), '-', 'Color',[0.75 0.2 0.2], 'DisplayName','\tau_x [N m]');
plot(tP, tau_onP(:,2), '-', 'Color',[0.2 0.75 0.2], 'DisplayName','\tau_y [N m]');
plot(tP, tau_onP(:,3), '-', 'Color',[0.2 0.2 0.75], 'DisplayName','\tau_z [N m]');
ylabel('[N m]');
yyaxis right
plot(tP, Tddot_onP,  '-','Color',[0.1 0.1 0.1], 'DisplayName','T̈_{cmd} [N/s^2]');
ylabel('[N/s^2]');
grid on; xlabel('t [s]'); title('Torques & T̈_{cmd} (ON)');
legend('Location','best'); if opts.MarkImpulses, mark_impulses(t0s); end
% bottom tile
nexttile; hold on;
plot(tP, rad2deg(etaacc_onP), '-', 'Color',[0.85 0.4 0.1], 'DisplayName','\etä_{cmd} [deg/s^2]');
plot(tP, rad2deg(theacc_onP), '-', 'Color',[0.1 0.6 0.85], 'DisplayName','\thetä_g_{cmd} [deg/s^2]');
grid on; xlabel('t [s]'); ylabel('[deg/s^2]'); title('Gimbal acceleration commands (ON)');
legend('Location','best'); if opts.MarkImpulses, mark_impulses(t0s); end
title(tlo2,'u_{QC} (M\_block \rightarrow QC dynamics) — ON');
saveplot(gcf, outDir, "u_QC_timeline__on", opts); close(f);

% ---------- return ----------
S = struct();
S.outDir = outDir;  
S.seqToken = seqToken;  
S.evalText = evalText;
S.metrics = struct('MAE3_off',MAE3_off,'MAE3_on',MAE3_on, ...
    'MAE_att_off',MAE_att_off,'MAE_att_on',MAE_att_on,'MAE_pov_off',MAE_pov_off,'MAE_pov_on',MAE_pov_on, ...
    'improvements',struct('pos',imp_pos,'alt',imp_alt,'attR',imp_attR,'attP',imp_attP,'attY',imp_attY,'pov',imp_pov));

% ---------- also save state histories for animation ----------
S.states = struct();
S.states.t     = t;
S.states.Xfw   = Xfw;
S.states.XqcOn = qc_on.X;
S.states.XqcOff= qc_off.X;

end  % ======================= end main =======================

% ======================= SUBFUNCTIONS (non-nested) =======================
function [params, initial_fw_state, initial_qc_state] = local_load_params()
    % Isolate param.m so any 'clearvars' inside it can't clear our workspace.
    param;  % must define: params, initial_fw_state, initial_qc_state
    if ~exist('params','var') || ~exist('initial_fw_state','var') || ~exist('initial_qc_state','var')
        error('param.m did not define required variables.');
    end
end

function S = run_qc_once(mode, params, x0_qc, ulog, t, dt, flags)
    N = numel(t);
    X = zeros(N,19); X(1,:) = x0_qc(:).';
    % logs
    Tlog=nan(N,1); nuTlog=nan(N,1); Tstar=nan(N,1); Tdd=nan(N,1);
    tau=nan(N,3); rcondBtB=nan(N,1); b3z=nan(N,1); lam=nan(N,1);
    wpsi=nan(N,1); gamma=nan(N,1); g1=nan(N,1); g2=nan(N,1);
    theta_err=nan(N,1); cos_gamma=nan(N,1);
    etaacc_log=nan(N,1); theacc_log=nan(N,1);

    for k=1:N-1
        xq = X(k,:).'; u  = ulog(k,:).';
        if k>1, params.resetFW = false; else, params.resetFW = true; params.fw_init_state = evalin('caller','initial_fw_state'); end
        [Tdd_cmd, trx, try_, trz, etaacc_cmd, theacc_cmd, ghost] = M_block(u, xq, params, mode); %#ok<ASGLU>
        % diagnostics
        if isfield(ghost,'diag')
            D = ghost.diag;
            if isfield(D,'T_star'),    Tstar(k)=D.T_star; end
            if isfield(D,'rcondBtB'),  rcondBtB(k)=D.rcondBtB; end
            if isfield(D,'b3z'),       b3z(k)=D.b3z; end
            if isfield(D,'lambda'),    lam(k)=D.lambda; end
            if isfield(D,'wpsi'),      wpsi(k)=D.wpsi; end
            if     isfield(D,'gamma'),       gamma(k)=D.gamma;
            elseif isfield(D,'gamma_fb'),    gamma(k)=D.gamma_fb;
            elseif isfield(D,'gamma_blend'), gamma(k)=D.gamma_blend;
            end
            if isfield(D,'g1'), g1(k)=D.g1; end
            if isfield(D,'g2'), g2(k)=D.g2; end
            if isfield(D,'theta_err'),  theta_err(k)=D.theta_err; end
            if isfield(D,'cos_gamma'),  cos_gamma(k)=D.cos_gamma; end
        end
        % logs to dynamics
        Tdd(k)=Tdd_cmd;  tau(k,:)=[trx, try_, trz];
        etaacc_log(k)=etaacc_cmd; theacc_log(k)=theacc_cmd;

        % integrate QC dynamics
        k1 = quadGimbalDynamics(xq,                Tdd_cmd,trx,try_,trz,etaacc_cmd,theacc_cmd,params);
        k2 = quadGimbalDynamics(xq+0.5*dt*k1,     Tdd_cmd,trx,try_,trz,etaacc_cmd,theacc_cmd,params);
        k3 = quadGimbalDynamics(xq+0.5*dt*k2,     Tdd_cmd,trx,try_,trz,etaacc_cmd,theacc_cmd,params);
        k4 = quadGimbalDynamics(xq+dt*k3,         Tdd_cmd,trx,try_,trz,etaacc_cmd,theacc_cmd,params);
        xqn = xq + (dt/6)*(k1+2*k2+2*k3+k4);
        xqn(14) = min(max(xqn(14), params.qc_T_min), params.qc_T_max);
        if flags.watchdog && ~all(isfinite(xqn))
            warning('[QC %s] state non-finite @k=%d', upper(mode),k);
            if flags.stop_on_bad, error('QC bad state'); end
            xqn(~isfinite(xqn)) = 0;
        end
        qq = xqn(7:10); nq = norm(qq); qq = qq/max(nq,1e-12); xqn(7:10)=qq;
        X(k+1,:) = xqn.';
        Tlog(k)=X(k,14); nuTlog(k)=X(k,15);
    end
    % last samples continuity
    Tlog(end)=X(end,14); nuTlog(end)=X(end,15);
    Tstar(end)=Tstar(end-1); Tdd(end)=Tdd(end-1);
    rcondBtB(end)=rcondBtB(end-1); b3z(end)=b3z(end-1); lam(end)=lam(end-1);
    wpsi(end)=wpsi(end-1); gamma(end)=gamma(end-1);
    g1(end)=g1(end-1); g2(end)=g2(end-1);
    theta_err(end)=theta_err(end-1); cos_gamma(end)=cos_gamma(end-1);
    etaacc_log(end)=etaacc_log(end-1); theacc_log(end)=theacc_log(end-1);

    S = struct('X',X,'T',Tlog,'nuT',nuTlog,'T_star',Tstar,'Tddot',Tdd, ...
               'tau',tau,'rcondBtB',rcondBtB,'b3z',b3z,'lambda',lam, ...
               'wpsi',wpsi,'gamma',gamma,'g1',g1,'g2',g2, ...
               'theta_err',theta_err,'cos_gamma',cos_gamma, ...
               'etaacc',etaacc_log,'theacc',theacc_log);
end

function [phi_qc,theta_qc,psi_qc,phi_cam,theta_cam,psi_cam,pov_err_deg] = derive_qc_views(Xqc, Xfw, ~)
    N = size(Xqc,1);
    phi_qc=zeros(N,1); theta_qc=zeros(N,1); psi_qc=zeros(N,1);
    phi_cam=zeros(N,1);theta_cam=zeros(N,1);psi_cam=zeros(N,1); pov_err_deg=zeros(N,1);
    Rx = @(a)[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
    Ry = @(a)[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
    for k=1:N
        qBN_q = Xqc(k,7:10).'; R_BN_q = quatNED_to_body(qBN_q); R_NB_q = R_BN_q.';
        [phi_qc(k), theta_qc(k), psi_qc(k)] = euler321_from_RNB(R_NB_q);
        qBN_f = Xfw(k,7:10).'; R_BN_f = quatNED_to_body(qBN_f); R_NB_f = R_BN_f.';
        eta = Xqc(k,16);  thg = Xqc(k,17);
        R_CB = Rx(eta)*Ry(thg);
        R_CN = R_CB.' * R_BN_q; R_NC = R_CN.';
        [phi_cam(k), theta_cam(k), psi_cam(k)] = euler321_from_RNB(R_NC);
        pov_err_deg(k) = rad2deg( rot_err(R_CN, R_BN_f) );
    end
end

function [phi,theta,psi] = euler321_from_RNB(R)
    theta = -asin( clamp(R(3,1), -1, 1) );
    phi   = atan2( R(3,2), R(3,3) );
    psi   = atan2( R(2,1), R(1,1) );
end

function y = clamp(x,a,b), y = min(max(x,a), b); end
function mark_impulses(t0s)
    % Draw impulse markers only if explicitly requested in opts.
    % Make them unobtrusive and *behind* plots on releases that support layering.
    for ii=1:numel(t0s)
        xl = xline(t0s(ii),'Color',[0.65 0.65 0.65], 'LineStyle',':', 'Alpha',0.15, ...
                   'HandleVisibility','off');                                      % ConstantLine
        if isprop(xl,'Layer'), xl.Layer = "bottom"; end
    end
end
function R = quatNED_to_body(q)
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);
    R = [ 1-2*(q2^2+q3^2),   2*(q1*q2 + q0*q3),  2*(q1*q3 - q0*q2);
          2*(q1*q2 - q0*q3), 1-2*(q1^2+q3^2),    2*(q2*q3 + q0*q1);
          2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1),  1-2*(q1^2+q2^2) ];
end
function ang = rot_err(R1,R2), s = (trace(R1.'*R2)-1)/2; s = min(max(s,-1),1); ang = acos(s); end

function plotPOVaxis(tP, y_off, y_on, y_fw, ylab, titleStr, outDir, C_FW,C_OFF,C_ON, VARS, t0s, opts)
    for v = VARS
        f = figure('Name',titleStr + " — " + v); hold on;
        if v=="both"
            plot(tP, y_off,'--','Color',C_OFF);
            plot(tP, y_on , '-.','Color',C_ON);
            plot(tP, y_fw , '-','Color',C_FW);
            legend('QC off','QC on','FW','Location','best');
        else
            plot(tP, y_on , '-','Color',C_ON);
            plot(tP, y_fw , '-','Color',C_FW);
            legend('QC on','FW','Location','best');
        end
        grid on; xlabel('t [s]'); ylabel(ylab); title(titleStr);
        if opts.MarkImpulses, mark_impulses(t0s); end
        saveplot(gcf, outDir, sanitize(titleStr) + "__" + v, opts); close(gcf);
    end
end

function plotNEDaxis(tP, y_off, y_on, y_fw, ylab, titleStr, outDir, C_FW,C_OFF,C_ON, VARS, t0s, opts)
    for v = VARS
        f = figure('Name',titleStr + " — " + v); hold on;
        if v=="both"
            plot(tP, y_off,'--','Color',C_OFF);
            plot(tP, y_on , '-.','Color',C_ON);
            plot(tP, y_fw , '-','Color',C_FW);
            legend('QC off','QC on','FW','Location','best');
        else
            plot(tP, y_on , '-','Color',C_ON);
            plot(tP, y_fw , '-','Color',C_FW);
            legend('QC on','FW','Location','best');
        end
        grid on; xlabel('t [s]'); ylabel(ylab); title(titleStr);
        if opts.MarkImpulses, mark_impulses(t0s); end
        saveplot(gcf, outDir, sanitize(titleStr) + "__" + v, opts); close(gcf);
    end
end

function plotPOV_offFW(tP, y_off, y_fw, ylab, titleStr, outDir, C_FW, C_OFF, t0s, opts)
    f = figure('Name', titleStr); hold on;
    plot(tP, y_off,'--','Color',C_OFF,'DisplayName','QC off');
    plot(tP, y_fw , '-','Color',C_FW ,'DisplayName','FW');
    grid on; xlabel('t [s]'); ylabel(ylab); title(titleStr);
    legend('Location','best');
    if opts.MarkImpulses, mark_impulses(t0s); end
    saveplot(gcf, outDir, sanitize(titleStr) + "__off_vs_FW", opts);
    close(f);
end

function saveplot(figH, outDir, stem, opts)
    % keep grid subtle; ensure visibility ordering; move legends off data
    axs = findall(figH,'Type','axes','-not','Tag','legend');
    for ax = reshape(axs,1,[])
        grid(ax,'on');   % idempotent; ensures grid even if a block forgot
        bringNonSolidToFront(ax);
        % If any ConstantLine (xline/yline) exists, push it behind data (R2024a+)
        cl = findall(ax,'Type','ConstantLine');
        for h = reshape(cl,1,[])
            if isprop(h,'Layer'), h.Layer = "bottom"; end
        end
    end
    if opts.LegendsOutside
        hL = findall(figH,'Type','Legend');
        for lgd = reshape(hL,1,[])
            lgd.Location = 'bestoutside';  % official option to avoid overlap
            lgd.Box = 'off';
        end
    end
    if opts.SavePNG
        exportgraphics(figH, fullfile(outDir, char(stem + ".png")), 'Resolution', 300);
    end
    if opts.SavePDF
        exportgraphics(figH, fullfile(outDir, char(stem + ".pdf")));
    end
end

function bringNonSolidToFront(ax)
    % Draw dashed/dotted lines on top of solid lines so overlaps stay visible.
    ch = ax.Children;
    if isempty(ch), return; end
    isLine = arrayfun(@(h) isa(h,'matlab.graphics.chart.primitive.Line'), ch);
    lines  = ch(isLine); if isempty(lines), return; end
    ls = get(lines,'LineStyle'); if ~iscell(ls), ls = {ls}; end
    solid    = lines(strcmp(ls,'-'));
    nonSolid = lines(~strcmp(ls,'-'));
    ax.Children = [ch(~isLine); solid; nonSolid];   % permutation of existing children only
end

function tok = seq_to_token(seq)
    n = size(seq,1); s = strings(n,1);
    for i=1:n
        ch = seq(i,1); k = seq(i,2);
        signChar = "p"; if k < 0, signChar = "n"; end
        s(i) = string(ch) + signChar + string(abs(k));
    end
    tok = "seq_" + strjoin(s,"_");
end

function s = sanitize(nameStr)
    s = string(nameStr);
    s = replace(s, [" ","/","\",":","*","?","""","<",">","|","—","–"], "_");
    s = regexprep(s, "_+", "_");
end

function out = ifelse(tf,a,b), if tf, out=a; else, out=b; end
end

function v = get_or(S, field, mask, defaultVal)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        a = S.(field); if isvector(a), a = a(:); end
        v = a(mask, :);
    else
        v = defaultVal;
    end
end

function xyz = MAE_xyz(XqcP, XfwP)
    dx = XqcP(:,1)-XfwP(:,1); dy = XqcP(:,2)-XfwP(:,2); dz = XqcP(:,3)-XfwP(:,3);
    xyz = [mean(abs(dx)) mean(abs(dy)) mean(abs(dz))];
end

function write_text(fname, txt)
    fid = fopen(fname,'wt'); if fid<0, warning('Could not open %s for writing.', fname); return; end
    fprintf(fid,'%s', txt); fclose(fid);
end

function txt = compose_eval_text(T_to, MAE3_off, MAE_xyz_off, MAX3_off, ...
                                      MAE_att_off,MAX_att_off,MAE_pov_off,MAX_pov_off, ...
                                      MAE3_on , MAE_xyz_on , MAX3_on , ...
                                      MAE_att_on , MAX_att_on , MAE_pov_on , MAX_pov_on , ...
                                      imp_pos,imp_alt,imp_attR,imp_attP,imp_attY,imp_pov, seqToken)
    b = string(newline);
    txt = "";
    txt = txt + sprintf('[Run] %s\n', seqToken);
    txt = txt + sprintf('[Window] 0..%.1f s\n\n', T_to);

    txt = txt + sprintf('[QC OFF vs FW]\n');
    txt = txt + sprintf('  3D pos MAE: %.3f m (x=%.3f y=%.3f z=%.3f); MAX=%.3f m\n', ...
        MAE3_off, MAE_xyz_off(1),MAE_xyz_off(2),MAE_xyz_off(3), MAX3_off);
    txt = txt + sprintf('  Att MAE [deg]: roll=%.2f (max=%.2f) pitch=%.2f (max=%.2f) yaw=%.2f (max=%.2f)\n', ...
        MAE_att_off(1),MAX_att_off(1),MAE_att_off(2),MAX_att_off(2),MAE_att_off(3),MAX_att_off(3));
    txt = txt + sprintf('  POV MAE [deg]: %.2f (max=%.2f)\n\n', MAE_pov_off, MAX_pov_off);

    txt = txt + sprintf('[QC ON  vs FW]\n');
    txt = txt + sprintf('  3D pos MAE: %.3f m (x=%.3f y=%.3f z=%.3f); MAX=%.3f m\n', ...
        MAE3_on, MAE_xyz_on(1),MAE_xyz_on(2),MAE_xyz_on(3), MAX3_on);
    txt = txt + sprintf('  Att MAE [deg]: roll=%.2f (max=%.2f) pitch=%.2f (max=%.2f) yaw=%.2f (max=%.2f)\n', ...
        MAE_att_on(1),MAX_att_on(1),MAE_att_on(2),MAX_att_on(2),MAE_att_on(3),MAX_att_on(3));
    txt = txt + sprintf('  POV MAE [deg]: %.2f (max=%.2f)\n\n', MAE_pov_on, MAX_pov_on);

    txt = txt + sprintf('[IMPROVEMENT ON vs OFF] (positive is better)\n');
    txt = txt + sprintf('  3D pos MAE:   %+5.2f %%\n', imp_pos);
    txt = txt + sprintf('  Alt MAE:      %+5.2f %%\n', imp_alt);
    txt = txt + sprintf('  Att MAE φ/θ/ψ:%+6.2f / %+6.2f / %+6.2f %%\n', imp_attR,imp_attP,imp_attY);
    txt = txt + sprintf('  POV MAE:      %+5.2f %%\n', imp_pov);
    txt = char(txt + b);
end
