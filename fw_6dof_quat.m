function xdot = fw_6dof_quat(state, thrust, elevator, aileron, rudder, params)
% FW_6DOF_QUAT  Rigid-body 6-DoF fixed-wing with wind-axis aero
%
% State (13x1): [x y z u v w q0 q1 q2 q3 p q r]'
%   x,y,z   : position in NED (m)  (Z down)
%   u,v,w   : body translational velocities (m/s)  (inertial-relative)
%   q0..q3  : quaternion q_BN (N->Body), scalar-first (Hamilton)
%   p,q,r   : body rates (rad/s)
%
% Inputs:
%   thrust   : propulsive thrust along +x_B (N)
%   elevator : elevator deflection (rad)
%   aileron  : aileron deflection  (rad)
%   rudder   : rudder deflection   (rad)
%
% Params (required fields; uses safe defaults if missing):
%   rho,S,b,c,g, wind_ned(3x1)
%   Ix,Iy,Iz,Ixz (preferred) OR J (3x3)
%   Aero coeffs: CL0, CL_alpha, CL_q, CL_de,
%                CD0, k (=induced), CDa, CD_q, CD_de,
%                CY_beta, CY_p, CY_r, CY_da, CY_dr,
%                Cl_beta, Cl_p, Cl_r, Cl_da, Cl_dr,
%                Cm0, Cm_alpha, Cm_q, Cm_de,
%                Cn_beta, Cn_p, Cn_r, Cn_da, Cn_dr
%   r_T (3x1): thrust offset from CG in body (m), optional (default [0;0;0])

% -------- unpack state
x = state(1); y = state(2); z = state(3); %#ok<NASGU>
u = state(4); v = state(5); w = state(6);
qBN = state(7:10);
p = state(11); q = state(12); r = state(13);

% -------- params & defaults
rho  = getdef(params,'rho',1.225);
S    = getdef(params,'S',0.25);
b    = getdef(params,'b',1.2);
c    = getdef(params,'c',S/max(b,eps));
g    = getdef(params,'g',9.81);
wind = getdef(params,'wind_ned',[0;0;0]);

% Inertia (allow Ixz coupling)
if all(isfield(params, {'Ix','Iy','Iz'}))
    Ix  = params.Ix; Iy = params.Iy; Iz = params.Iz;
    Ixz = getdef(params,'Ixz',0);
    J   = [Ix, 0, -Ixz;
           0,  Iy, 0;
          -Ixz,0,  Iz];
else
    J = getdef(params,'J',diag([0.045 0.065 0.085]));
end

% Aero coefficients (provide zeros if missing)
CL0   = getdef(params,'CL0',0);
CL_a  = getdef(params,'CL_alpha',0);
CL_q  = getdef(params,'CL_q',0);
CL_de = getdef(params,'CL_de',0);

CD0   = getdef(params,'CD0',0.01);
k_ind = getdef(params,'k',0.03);
CD_a  = getdef(params,'CDa',0);
CD_q  = getdef(params,'CD_q',0);
CD_de = getdef(params,'CD_de',0);

CY_b  = getdef(params,'CY_beta',0);
CY_p  = getdef(params,'CY_p',0);
CY_r  = getdef(params,'CY_r',0);
CY_da = getdef(params,'CY_da',0);
CY_dr = getdef(params,'CY_dr',0);

Cl_b  = getdef(params,'Cl_beta',0);
Cl_p  = getdef(params,'Cl_p',0);
Cl_r  = getdef(params,'Cl_r',0);
Cl_da = getdef(params,'Cl_da',0);
Cl_dr = getdef(params,'Cl_dr',0);

Cm0   = getdef(params,'Cm0',0);
Cm_a  = getdef(params,'Cm_alpha',0);
Cm_q  = getdef(params,'Cm_q',0);
Cm_de = getdef(params,'Cm_de',0);

Cn_b  = getdef(params,'Cn_beta',0);
Cn_p  = getdef(params,'Cn_p',0);
Cn_r  = getdef(params,'Cn_r',0);
Cn_da = getdef(params,'Cn_da',0);
Cn_dr = getdef(params,'Cn_dr',0);

r_T = getdef(params,'r_T',[0;0;0]);

% -------- rotations
R_BN = quatNED_to_body(qBN);  % N->B
R_NB = R_BN.';                % B->N

% -------- kinematics: position time-derivative in NED
v_b   = [u; v; w];
v_ned = R_NB * v_b;
x_dot = v_ned(1);
y_dot = v_ned(2);
z_dot = v_ned(3);

% -------- air-relative velocity in body
v_air_b = v_b - R_BN*wind;    % (R_BN*wind) is wind expressed in body
ua = v_air_b(1); va = v_air_b(2); wa = v_air_b(3);
Va = max(1e-3, norm(v_air_b));

alpha = atan2(wa, ua);                     % angle of attack
beta  = asin( clamp(va/Va, -1, 1) );       % sideslip

qbar = 0.5 * rho * Va^2;
p_hat = (b/(2*Va)) * p;
q_hat = (c/(2*Va)) * q;
r_hat = (b/(2*Va)) * r;

% -------- aerodynamic coefficients
CL = CL0 + CL_a*alpha + CL_q*q_hat + CL_de*elevator;
CD = CD0 + k_ind*CL^2 + CD_a*alpha + CD_q*q_hat + CD_de*elevator;
CY = CY_b*beta + CY_p*p_hat + CY_r*r_hat + CY_da*aileron + CY_dr*rudder;

Cl = Cl_b*beta + Cl_p*p_hat + Cl_r*r_hat + Cl_da*aileron + Cl_dr*rudder;
Cm = Cm0 + Cm_a*alpha + Cm_q*q_hat + Cm_de*elevator;
Cn = Cn_b*beta + Cn_p*p_hat + Cn_r*r_hat + Cn_da*aileron + Cn_dr*rudder;

% -------- wind-axis → body conversion for forces
% In wind axes, F_w = [ -D,  Y,  -L ] * qbar*S
F_w = qbar*S * [-CD; CY; -CL];

% Rotation from body->wind (C_b^w), then use its transpose for wind->body
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta);  sb = sin(beta);
C_b_w = [  ca*cb,  sb,    sa*cb;
          -ca*sb,  cb,   -sa*sb;
          -sa,     0,      ca  ];
R_wb = C_b_w.';               % wind -> body
F_aero_b = R_wb * F_w;

% -------- aerodynamic moments in body
M_aero_b = qbar*S * [ b*Cl; c*Cm; b*Cn ];

% -------- thrust force & moment
F_thrust_b = [thrust; 0; 0];
M_thrust_b = cross(r_T(:), F_thrust_b);

% -------- gravity in body
F_grav_b = R_BN * [0;0; g*getdef(params,'m',1.0)]; % using FW mass if provided
m_fw = getdef(params,'m',1.0);

% -------- total body forces & accelerations
F_b = F_aero_b + F_thrust_b + F_grav_b;

% body translational dynamics:  v̇_b = (1/m)F_b - ω×v_b
omega = [p;q;r];
v_dot_b = (1/m_fw) * F_b - cross(omega, v_b);

u_dot = v_dot_b(1);
v_dot = v_dot_b(2);
w_dot = v_dot_b(3);

% -------- attitude kinematics (quaternion)
Omega = [ 0   -p   -q   -r;
          p    0    r   -q;
          q   -r    0    p;
          r    q   -p    0 ];
q_dot = 0.5 * Omega * qBN;     % (normalize after integration in driver)

% -------- rotational dynamics
M_b = M_aero_b + M_thrust_b;
omega_dot = J \ ( M_b - cross(omega, J*omega) );
p_dot = omega_dot(1);
q_dotb = omega_dot(2);
r_dot = omega_dot(3);

% -------- assemble derivative
xdot = zeros(13,1);
xdot(1)  = x_dot;
xdot(2)  = y_dot;
xdot(3)  = z_dot;
xdot(4)  = u_dot;
xdot(5)  = v_dot;
xdot(6)  = w_dot;
xdot(7:10) = q_dot;
xdot(11) = p_dot;
xdot(12) = q_dotb;
xdot(13) = r_dot;

end

% ===================== helpers ==========================
function R = quatNED_to_body(q)
% q = [q0 q1 q2 q3] (scalar-first), mapping NED -> Body (Hamilton)
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
R = [ 1-2*(q2^2+q3^2),   2*(q1*q2 + q0*q3),  2*(q1*q3 - q0*q2);
      2*(q1*q2 - q0*q3), 1-2*(q1^2+q3^2),    2*(q2*q3 + q0*q1);
      2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1),  1-2*(q1^2+q2^2) ];
end

function y = clamp(x, a, b)
y = min(max(x,a), b);
end

function out = getdef(S, field, defaultVal)
if isfield(S, field), out = S.(field); else, out = defaultVal; end
end
