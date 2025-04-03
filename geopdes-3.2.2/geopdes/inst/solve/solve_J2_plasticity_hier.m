% SOLVE_LINEAR_ELASTICITY: Solve an elasto-plastic problem on a NURBS domain.
%                          For a planar domain it is the plane strain model.
%
% The function solves the elasto-plastic problem using von Mises criterion
%
%      - div (sigma(u)) = f    in Omega = F((0,1)^n)
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% with   sigma(u) = mu*(grad(u) + grad(u)^t) + lambda*div(u)*I.
%
%   u:          displacement vector
%   sigma:      Cauchy stress tensor
%   lambda, mu: Lame' parameters
%   I:          identity tensor
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_J2_plasticity (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - press_sides:  sides with pressure boundary condition (may be empty)
%    - symm_sides:   sides with symmetry boundary condition (may be empty)
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - nload:      number of load steps increment
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete basis functions (see sp_vector)
%  u:        the computed degrees of freedom
%
% See also EX_LIN_ELAST_HORSESHOE for an example.
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [u, eps_pl, sigma] = solve_J2_plasticity_hier (problem_data, method_data, adaptivity_data, hspace, hmsh, load_multiplier, u, eps_pl, sigma)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end
data_names = fieldnames (adaptivity_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= adaptivity_data.(data_names{iopt});']);
end

% Assemble the matrices
rhs_shape    = op_f_v_hier (hspace, hmsh, f);

% Apply Neumann boundary conditions
for iside = nmnn_sides
    % Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
    gside = @(varargin) g(varargin{:},iside);
    dofs = hspace.boundary(iside).dofs;
    rhs_shape(dofs) = rhs_shape(dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), gside);
end

% Apply pressure conditions
for iside = problem_data.press_sides
    pside = @(varargin) problem_data.p(varargin{:},iside);
    dofs = hspace.boundary(iside).dofs;
    rhs_shape(dofs) = rhs_shape(dofs) - op_pn_v_hier (hspace.boundary(iside), hmsh.boundary(iside), pside);
end

% Apply symmetry conditions
symm_dofs = [];
for iside = symm_sides
    msh_side = hmsh_eval_boundary_side (hmsh, iside);
    normal_comp = zeros (msh_side.rdim, msh_side.nqn * msh_side.nel);
    for idim = 1:msh_side.rdim
        normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
    end

    parallel_to_axes = false;
    for ind = 1:msh_side.rdim
        ind2 = setdiff (1:msh_side.rdim, ind);
        if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
            symm_dofs = union (symm_dofs, hspace.boundary(iside).dofs(hspace.boundary(iside).comp_dofs{ind}));
            parallel_to_axes = true;
            break
        end
    end
    if (~parallel_to_axes)
        error ('adaptivity_solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
    end
end

% Apply Dirichlet boundary conditions
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (hspace, hmsh, h, drchlt_sides);
int_dofs = setdiff (1:hspace.ndof, union (drchlt_dofs, symm_dofs));

%% Solve linear system
rhs = rhs_shape./nload .* load_multiplier;
iter = 0;
u(drchlt_dofs) = u_drchlt;

% Assemble tangent matrix
[K, internal_energy, eps_pl_new, sigma_new] = op_plsu_ev_hier (hspace, hspace, hmsh, u, eps_pl, sigma, mu_lame, kappa_lame, yield_stress);
external_energy = rhs(int_dofs) - K(int_dofs, drchlt_dofs) * u_drchlt;
res = internal_energy(int_dofs) - external_energy;
res_norm_0 = norm(res);
res_norm =res_norm_0;

% Newton-Raphson while loop
while res_norm/res_norm_0 > method_data.newton_tol && res_norm > method_data.newton_tol_abs && iter < method_data.newton_iter_max 

    % Solve the nonlinear system
    u_inc = - K(int_dofs, int_dofs) \ res;

    % Update solution vector
    u(int_dofs) = u(int_dofs) + u_inc;

    % Evaluate residuum
    [K, internal_energy, eps_pl_new, sigma_new] = op_plsu_ev_hier (hspace, hspace, hmsh, u, eps_pl, sigma, mu_lame, kappa_lame, yield_stress);
    external_energy = rhs(int_dofs) - K(int_dofs, drchlt_dofs) * u_drchlt;
    res = internal_energy(int_dofs) - external_energy;
    res_norm = norm(res)
    iter = iter+1

end % end N-R while loop
eps_pl = eps_pl_new;
sigma = sigma_new;


