%% ===========================  param.m  ===============================
% Quaternion-only setup for:
%   • fixed-wing dynamics   (fw_6dof_quat)            [state: 13]
%   • quad + 2-axis gimbal  (quadGimbalDynamics)      [state: 19]
%   • mapping / M_block (NDI/DFL with a ghost FW)
%
% Frames: NED (Z down). Quaternions scalar-first [q0 q1 q2 q3].
% Convention for both plants: q_BN (rotation NED->Body).
% ----------------------------------------------------------------------

clearvars; clc;

%% ───────────────────────── Universal constants ───────────────────────
params.g        = 9.81;                 % +Z_ned is down
params.rho      = 1.225;                % kg/m^3
params.dt       = 0.001;                % base step (s)
params.wind_ned = [0;0;0];              % m/s (override per scenario)

%% ───────────── Quad-copter + gimbal physical parameters ─────────────
% ~1 kg quad typical numbers (450-class)
params.qc_m    = 1.00;                                  % kg
params.qc_J    = diag([0.020, 0.020, 0.035]);           % kg·m^2 (Ixx,Iyy,Izz)
params.qc_l    = 0.20;                                  % m  (arm length)

% Actuator limits (QC plant)
params.qc_T_max   = 30;                                 % N (total thrust max)
%params.qc_T_min   = 0;       % N
params.dfl.T_floor_frac = 0.10;                    % 10% of mg (was 60%)
params.qc_T_min         = params.dfl.T_floor_frac * params.qc_m * params.g;

params.qc_tau_max = [1.2; 1.2; 1.2];                    % N·m [roll; pitch; yaw]

% Simple quadratic body drag for QC (per-axis Cd*A, m^2) — keep 0 to match paper
params.qc_CdA = [0, 0, 0];

% ---- Aliases used by M_block (don’t touch FW fields below) ----
params.m_qc     = params.qc_m;      % QC mass for M_block
params.J_qc     = params.qc_J;      % QC inertia for M_block
params.tau_max  = params.qc_tau_max;% torque saturation used by M_block

% Thrust-channel 2nd-order tracker (T, nu_T := Tdot)
params.kT0 = 120;     % "stiffness"
params.kT1 =  30;     % "damping"

% --- Safety attitude feedback (QC -> ghost) gains for M_block ---
params.fb.kR  = [0.25; 0.25; 0.35];   % N·m per rad  (roll, pitch, yaw)
params.fb.kOm = [0.03; 0.03; 0.04];   % N·m per rad/s
params.fb.theta_act = deg2rad(10);    % [rad] start blending around 10 deg att error
params.fb.max_blend = 0.7;            % 0..1, max fraction of PD you allow


%% ───────────────────────── 2-Axis gimbal settings ────────────────────
% Kinematics: R_CB = R_x(eta) * R_y(theta_g)
params.gimbal_axes    = 2;                             % roll–pitch (no yaw)
params.gim_rate_max   = deg2rad([360; 360]);           % very high (near-ideal)
params.gim_ang_max    = deg2rad([60; 45]);             % [roll; pitch] rad (±)
params.kg0_roll       = 80.0;                          % PD gains → joint acc cmd
params.kg1_roll       = 20.0;
params.kg0_pitch      = 80.0;
params.kg1_pitch      = 20.0;

% Ideal gimbal option (snap to command, no dynamics/limits):
params.gimbal_ideal   = false;   % set true to make it "perfect/fixed to angle"

%% ───────────── Fixed-wing (FW) inertial & aero parameters ────────────
% Names match fw_6dof_quat() expectations (FW plant ONLY).
params.m  = 1.00;                                       % kg (FW mass)
params.J  = diag([0.045, 0.065, 0.085]);                % kg·m^2 (Ixx,Iyy,Izz)
params.S  = 0.25;                                       % m^2
params.b  = 1.20;                                       % m
params.c  = params.S / params.b;                        % m

% Lift
params.CL0  = 0.40;
params.CLa  = 7.5;      % 1/rad
params.CLq  = 0.0;
params.CLde = 2.0;

% Drag
params.CD0  = 0.005;
params.CDa  = 0.0;
params.CDq  = 0.0;
params.CDde = 0.0;
params.CDk  = 0.03;     % induced term k*CL^2

% Side force
params.CYb  = -1.5;
params.CYp  = 0.0;
params.CYr  = 0.0;
params.CYda = 0.0;
params.CYdr = 0.0;

% Roll moment
params.Clb  = -0.12;
params.Clp  = -0.5;
params.Clr  = 0.25;
params.Clda = 0.012;
params.Cldr = 0.0;

% Pitch moment
params.Cm0  = 0.00;
params.Cma  = -0.50;
params.Cmq  = -10.0;
params.Cmde = -0.24;

% Yaw moment
params.Cnb  = 0.20;
params.Cnp  = -0.06;
params.Cnr  = -0.20;
params.Cnda = 0.00;
params.Cndr = 0.02;

%% ───────────── Fixed-wing compat map for fw_6dof_quat() ─────────────
% Inertia scalars expected by FW plant
params.Ix  = params.J(1,1);
params.Iy  = params.J(2,2);
params.Iz  = params.J(3,3);
params.Ixz = 0.0;                 % kg·m^2, small/symmetric UAV ≈ 0

% Induced drag factor used by FW plant (CD = CD0 + k*CL^2)
params.k   = params.CDk;

% Lift (names expected by FW plant)
params.CL_alpha = params.CLa;     % 1/rad
params.CL_q     = params.CLq;
params.CL_de    = params.CLde;

% Side force
params.CY_beta  = params.CYb;
params.CY_p     = params.CYp;
params.CY_r     = params.CYr;
params.CY_da    = params.CYda;
params.CY_dr    = params.CYdr;

% Roll moment
params.Cl_beta  = params.Clb;
params.Cl_p     = params.Clp;
params.Cl_r     = params.Clr;
params.Cl_da    = params.Clda;
params.Cl_dr    = params.Cldr;

% Pitch moment
params.Cm_alpha = params.Cma;
params.Cm_q     = params.Cmq;
params.Cm_de    = params.Cmde;

% Yaw moment
params.Cn_beta  = params.Cnb;
params.Cn_p     = params.Cnp;
params.Cn_r     = params.Cnr;
params.Cn_da    = params.Cnda;
params.Cn_dr    = params.Cndr;


% Optional: thrust-line offset (body, from CG) for FW thrust moment
if ~isfield(params,'r_T'); params.r_T = [0;0;0]; end

%Vertical guard
if ~isfield(params,'dfl'), params.dfl = struct(); end
params.dfl.T_star_vertical_guard = false;   % disable the vertical-guard

%% ───────────── FW open-loop helper (thrust limits, etc.) ─────────────
params.max_thrust = 25;    % N
params.min_thrust = 0.001;     % N

%% ───────────── Source for FW initial state / trace (optional) ────────
use_fw_trace = false;       % if true and file exists, take initial FW state from file

filename = 'FW20SecondsRun.xlsx';
if use_fw_trace && isfile(filename)
    tbl = readtable(filename,'Sheet',1);
    names = [{'time'}, arrayfun(@(k)sprintf('outputStates%d',k),1:13,'uni',false)];
    n = min(numel(names), width(tbl));
    tbl.Properties.VariableNames(1:n) = names(1:n);

    % [x y z u v w q0 q1 q2 q3 p q r]   (FW convention here is body->NED)
    S = tbl{:, 2:14};   % [N x 13]
    if size(S,2) ~= 13
        error('Expected 13 state columns (outputStates1..13). Found %d.', size(S,2));
    end

    % Quick quaternion norm check (cols 7:10)
    qn = sqrt(sum(S(:,7:10).^2,2));
    if any(abs(qn-1) > 1e-3)
        warning('FW quaternions deviate from unit norm by >1e-3 at %d samples.', sum(abs(qn-1) > 1e-3));
    end

    FW_States = [tbl.time, S];                         % [time, 13 states]
    initial_fw_state = FW_States(1, 2:14).';           % 13x1
else
    % No file → build a straight-and-level initial FW state at 15 m/s North.
    vN0 = 15; vE0 = 0; vD0 = 0;
    x0 = 0; y0 = 0; z0 = 50;                            % 50 m AGL, Z down
    q_BN0 = [1;0;0;0];                                  % N->B identity (body fwd = +N)
    R_NB0 = eye(3);
    uvw0  = R_NB0.' * [vN0; vE0; vD0];                  % body velocities
    pqr0  = [0;0;0];

    % FW state uses q_BI or q_BN? We standardize on q_BN everywhere.
    initial_fw_state = [x0; y0; z0; uvw0; q_BN0; pqr0]; % 13x1
    FW_States = [0, initial_fw_state.'];                % stub row
end

%% ──────────── Initial states (FW and QC+gimbal, quaternion) ──────────
% FW initial: [x y z u v w q0 q1 q2 q3 p q r]'
% (Our FW dynamics will read q as q_BN, consistent with QC.)
fw = initial_fw_state;
x0 = fw(1); y0 = fw(2); z0 = fw(3);
u0 = fw(4); v0 = fw(5); w0 = fw(6);
q_BN_fw = fw(7:10);                                     % N->Body, scalar-first
pqr0    = fw(11:13);

% Body->NED (for convenience) to compute v_ned0:
% Prefer quat2rotm if quat2dcm missing.
if exist('quat2dcm','file')
    R_NB0 = quat2dcm([q_BN_fw(1), -q_BN_fw(2), -q_BN_fw(3), -q_BN_fw(4)]); % conj gives B->N
else
    R_BN  = [ q_BN_fw(1), q_BN_fw(2), q_BN_fw(3), q_BN_fw(4) ];
    R_NB0 = quat2rotm([R_BN(1), -R_BN(2), -R_BN(3), -R_BN(4)]);
end
v_ned0  = R_NB0 * [u0; v0; w0];

% Thrust dynamic extension initial values for QC
m_qc = params.qc_m;
T0   = min(max(m_qc * params.g, params.qc_T_min), params.qc_T_max); % hover-ish
nu_T0 = 0;

% 2-axis gimbal angles & rates
eta0 = 0; theta_g0 = 0; etadot0 = 0; thetadot0 = 0;

% QC+g state vector (19x1):
% [x y z vN vE vD q0 q1 q2 q3 p q r T nu_T eta theta_g eta_dot theta_g_dot]
initial_qc_state = [ ...
    x0; y0; z0; ...
    v_ned0(:); ...
    q_BN_fw; ...
    pqr0; ...
    T0; nu_T0; ...
    eta0; theta_g0; etadot0; thetadot0 ];

%% ─────────── NDI/DFL controller gains (M_block) ──────────────────────
% 4th-order position-chain gains (per-axis allowed; here scalar → same for xyz)
params.k0_pos = 8.0;   % position
params.k1_pos = 8.0;   % velocity
params.k2_pos = 5.0;   % acceleration
params.k3_pos = 3.0;   % jerk

% Yaw 2nd-order chain gains
params.kpsi0  = 6.0;   % yaw position gain
params.kpsi1  = 3.0;   % yaw rate gain

% Ghost FW state for M_block (and optional reset flag)
params.fw_init_state = initial_fw_state;  % 13x1 (internal ghost copy)
params.resetFW       = false;             % set true to reinit ghost FW next call

%% ─────────── Export to workspace ─────────────────────────────────────
assignin('base','params',            params);
assignin('base','FW_States',         FW_States);
assignin('base','initial_fw_state',  initial_fw_state);   % 13×1
assignin('base','initial_qc_state',  initial_qc_state);   % 19×1

% (Optional) state index maps to avoid magic numbers in models
idx_fw = struct('x',1,'y',2,'z',3,'u',4,'v',5,'w',6,'q0',7,'q1',8,'q2',9,'q3',10,'p',11,'q',12,'r',13);
idx_qc = struct('x',1,'y',2,'z',3,'vN',4,'vE',5,'vD',6,'q0',7,'q1',8,'q2',9,'q3',10,'p',11,'q',12,'r',13, ...
                'T',14,'nu_T',15,'eta',16,'theta_g',17,'eta_dot',18,'theta_g_dot',19);
assignin('base','idx_fw', idx_fw);
assignin('base','idx_qc', idx_qc);

% Optional: create Simulink bus for params
if exist('Simulink.Bus.createObject','file')
    paramsBusInfo = Simulink.Bus.createObject(params);
    paramsBus     = evalin('base', paramsBusInfo.busName);
    assignin('base','paramsBus', paramsBus);
end

fprintf(['[param] Params loaded. FW_States [%dx14], ', ...
         'initial FW/QC states exported (FW 13x1, QC 19x1).\n'], ...
         size(FW_States,1));

% (Optional) load additional datasets
if exist('loadFixedWingData','file')==2
    loadFixedWingData;
end
