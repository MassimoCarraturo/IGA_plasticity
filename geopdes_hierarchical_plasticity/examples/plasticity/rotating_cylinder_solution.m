function [sigma_r, sigma_t, u_r, r_eval] = rotating_cylinder_solution(E, nu, rho_omega2, R_i, R_o, n_points)
% ROTATING_CYLINDER_SOLUTION  Elastic analytical solution for a thick-walled
% rotating cylinder in plane strain.
%
% The stress field for a spinning cylinder with traction-free inner and
% outer surfaces satisfies the Navier equation with centrifugal body force
%   f_r = rho * omega^2 * r.
%
% Plane-strain closed-form (Timoshenko & Goodier, §42):
%   sigma_r = A + B/r^2 - C3 * r^2
%   sigma_t = A - B/r^2 - C1 * r^2
%   u_r     = (1/E') * [ (1-nu')*A*r - (1+nu')*B/r - C_u*r^3 ]
%
% where:
%   C3 = (3 - 2*nu) / (8*(1 - nu)) * rho*omega^2
%   C1 = (1 + 2*nu) / (8*(1 - nu)) * rho*omega^2
%
% BCs: sigma_r(R_i) = 0,  sigma_r(R_o) = 0
%
% Input:
%   E, nu         - Young's modulus [MPa], Poisson's ratio
%   rho_omega2    - product rho*omega^2 [MPa/mm^2] (body force coefficient)
%   R_i, R_o      - inner and outer radius [mm]
%   n_points      - number of evaluation points along radius
%
% Output:
%   sigma_r  - radial stress [MPa],   size [1 x n_points]
%   sigma_t  - tangential stress [MPa], size [1 x n_points]
%   u_r      - radial displacement [mm], size [1 x n_points]
%   r_eval   - evaluation radii [mm],   size [1 x n_points]

r_eval = linspace(R_i, R_o, n_points);

% Plane-strain coefficients
C3 = (3 - 2*nu) / (8*(1 - nu)) * rho_omega2;
C1 = (1 + 2*nu) / (8*(1 - nu)) * rho_omega2;

% BCs: sigma_r(R_i) = 0 => A + B/R_i^2 - C3*R_i^2 = 0
%       sigma_r(R_o) = 0 => A + B/R_o^2 - C3*R_o^2 = 0
% Solve 2x2 system: [1, 1/R_i^2; 1, 1/R_o^2] * [A; B] = [C3*R_i^2; C3*R_o^2]
mat = [1, 1/R_i^2; 1, 1/R_o^2];
rhs_vec = [C3*R_i^2; C3*R_o^2];
AB = mat \ rhs_vec;
A = AB(1);
B = AB(2);

% Stress field
sigma_r = A + B ./ r_eval.^2 - C3 * r_eval.^2;
sigma_t = A - B ./ r_eval.^2 - C1 * r_eval.^2;

% Displacement (plane strain): use effective E' = E/(1-nu^2), nu' = nu/(1-nu)
% u_r = r/E * [ (1-nu)*sigma_t - nu*sigma_r ] + plane-strain correction
% More directly from the Lamé solution:
%   eps_r = (1/E)*[(1-nu^2)*sigma_r - nu*(1+nu)*sigma_t] + (rho_omega2 term)
% Or use the standard result:
%   u_r = (1+nu)/E * [ (1-2*nu)*A*r - B/r - (1-2*nu)*C3*r^3/(3-2*nu) ...
%                       + ?? ]
% Let me use the direct plane-strain formula:
%   eps_theta = u_r / r
%   eps_r     = du_r/dr
%   sigma_r = E/((1+nu)(1-2nu)) * [(1-nu)*eps_r + nu*eps_theta]
%   sigma_t = E/((1+nu)(1-2nu)) * [nu*eps_r + (1-nu)*eps_theta]
%
% From sigma_t = A - B/r^2 - C1*r^2 and sigma_r = A + B/r^2 - C3*r^2:
%   eps_theta = u_r/r
% Working backwards from the stress-displacement relations in plane strain:
%   u_r = r * [ (1+nu)/E * ( (1-nu)*sigma_t - nu*sigma_r ) ]
%       = r * (1+nu)/E * [ (1-nu)*(A - B/r^2 - C1*r^2) - nu*(A + B/r^2 - C3*r^2) ]
%       = (1+nu)/E * [ (1-2*nu)*A*r - B/r + r^3 * (-C1*(1-nu) + nu*C3) ]

coeff_r3 = -C1*(1-nu) + nu*C3;
u_r = (1+nu)/E * ( (1-2*nu)*A*r_eval - B./r_eval + coeff_r3 * r_eval.^3 );

end
