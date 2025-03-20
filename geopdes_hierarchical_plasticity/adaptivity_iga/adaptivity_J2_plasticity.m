% ADAPTIVITY_linear_elasticity: solve the linear elasticity problem with an
% adaptive isogeometric method based on hierarchical splines.
%
% [geometry, hmsh, hspace, u, solution_data] = adaptivity_linear_elasticity (problem_data, method_data, adaptivity_data, plot_data)
%
% XXXXXXXXXXXXX TEXT IS MISSING
%
% Copyright (C) 2017 Cesare Bracco, Rafael Vazquez
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

function [geometry, cell_hmsh, cell_hspace, cell_u, solution_data] = adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data, plot_data)

if (nargin == 3)
    plot_data = struct ('print_info', true, 'plot_hmesh', false, 'plot_discrete_sol', false);
end
if (~isfield (plot_data, 'print_info'))
    plot_data.print_info = true;
end
if (~isfield (plot_data, 'plot_hmesh'))
    plot_data.plot_hmesh = false;
end
if (~isfield (plot_data, 'plot_discrete_sol'))
    plot_data.plot_discrete_sol = false;
end

% Initialization of some auxiliary variables
if (plot_data.plot_hmesh)
    fig_mesh = figure;
end
if (plot_data.plot_discrete_sol)
    fig_sol = figure;
end
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+1;

% initialize space for solution
cell_hmsh = cell(method_data.nload,1);
cell_hspace = cell(method_data.nload,1);
cell_u = cell(method_data.nload,1);
cell_eps_pl = cell(method_data.nload,1); 

% Initialization of the hierarchical mesh and space
[hmsh, hspace, geometry] = adaptivity_initialize_vector (problem_data, method_data);

% store initial scalar space
[knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);
rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, geometry);
space_scalar =sp_bspline (knots, method_data.degree, msh);
clear knots zeta rule qn  qw msh
hspace_scalar   = hierarchical_space (hmsh, space_scalar, method_data.space_type, method_data.truncated, method_data.regularity);
% Initialization of problem variables
u = zeros (hspace.ndof,1);
eps_pl = cell(hmsh.nlevels,1);
for ilev = 1:hmsh.nlevels
    eps_pl{ilev} = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
end




%% Analysis

% LOAD INCREMENT LOOP
for iLoad = 1: method_data.nload

    % ADAPTIVE LOOP
    iter = 0;
    while (1)
        iter = iter + 1;

        if (plot_data.print_info)
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
        end
        if (~hspace_check_partition_of_unity (hspace, hmsh))
            disp('ERROR: The partition-of-the-unity property does not hold.')
            solution_data.flag = -1; break
        end

        % SOLVE 
        if (plot_data.print_info)
            disp('SOLVE:')
            fprintf('Number of elements: %d. Total DOFs: %d \n', hmsh.nel, hspace.ndof);
        end
        [u, eps_pl] = solve_J2_plasticity_hier (problem_data, method_data, adaptivity_data, hspace, hmsh, iLoad, u, eps_pl);

        % QI or L2 for eps_pl projection
        
        
        eps_pl_store = zeros(hspace_scalar.ndof, 6); % control variables
        eps_pl_store(:,:) = history_variable_projection_hier(hspace_scalar, hmsh, eps_pl,  method_data.type_projection);



        nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;

       

        % ESTIMATE
        if (plot_data.print_info); disp('ESTIMATE:'); end
        est = ones(hmsh.nel,1);% adaptivity_estimate_linear_el (u, hmsh, hspace, problem_data, adaptivity_data);
        est_elem{iter} = est.';
        gest(iter) = norm (est);
        if (plot_data.print_info); fprintf('Computed error estimate: %e \n', gest(iter)); end
        if (isfield (problem_data, 'graduex'))
            [err_h1(iter), err_l2(iter), err_h1s(iter),~,~,err_h1s_elem{iter}] = sp_h1_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex);
            if (plot_data.print_info); fprintf('Error in H1 seminorm = %e\n', err_h1s(iter)); end
        elseif (isfield (problem_data, 'uex'))
            err_l2(iter) = sp_l2_error (hspace, hmsh, u, problem_data.uex);
            if (plot_data.print_info); fprintf('Error in L2 norm = %e\n', err_l2(iter)); end
        end

        % STOPPING CRITERIA
        if (gest(iter) < adaptivity_data.tol)
            disp('Success: The error estimation reached the desired tolerance');
            solution_data.flag = 1; break
        elseif (iter == adaptivity_data.num_max_iter)
            disp('Warning: reached the maximum number of iterations')
            solution_data.flag = 2; break
        elseif (hmsh.nlevels >= adaptivity_data.max_level)
            disp('Warning: reached the maximum number of levels')
            solution_data.flag = 3; break
        elseif (hspace.ndof > adaptivity_data.max_ndof)
            disp('Warning: reached the maximum number of DOFs')
            solution_data.flag = 4; break
        elseif (hmsh.nel > adaptivity_data.max_nel)
            disp('Warning: reached the maximum number of elements')
            solution_data.flag = 5; break
        end

        % MARK
        if (plot_data.print_info); disp('MARK:'); end
        [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);
        if (plot_data.print_info)
            fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
            disp('REFINE:')
        end

        % REFINE
        hmsh_coarse = hmsh;
        hspace_coarse = hspace;
        hspace_scalar_coarse  = hspace_scalar;
        [hmsh, hspace, Cref] = adaptivity_refine (hmsh_coarse, hspace_coarse, marked, adaptivity_data);
        [~, hspace_scalar, ~] = adaptivity_refine (hmsh_coarse, hspace_scalar_coarse, marked, adaptivity_data);

        % refine variables
        u = Cref * u;
        eps_pl_store = Cref(1:2:end,1:2:end)*eps_pl_store; % control variables
        eps_pl = evaluate_at_quad_points(hmsh, hspace, eps_pl_store);


        % L2-projection
        % eps_pl = eps_pl_project_l2(hmsh_coarse, hspace_coarse, eps_pl, hmsh, hspace, Cref);
        % quasi-interpolant
        % eps_pl = eps_pl_quasi_interpolant(hmsh_coarse, hmsh, hspace, eps_pl);

        fprintf('\n');
    end % end adaptivity step

    solution_data.iter = iter;
    solution_data.gest = gest(1:iter);
    solution_data.ndof = ndof(1:iter);
    solution_data.nel  = nel(1:iter);
    if (exist ('err_h1s', 'var'))
        solution_data.err_h1s = err_h1s(1:iter);
        solution_data.err_h1 = err_h1(1:iter);
        solution_data.err_l2 = err_l2(1:iter);
        solution_data.err_h1s_elem=err_h1s_elem;
    end

    % store solution

    cell_hmsh{iLoad} = hmsh;
    cell_hspace{iLoad} = hspace;
    cell_u{iLoad} = u;
    cell_eps_pl{iLoad} = eps_pl_store;


end % end load step

end






%--------------------------------------------------------------------------
% subroutines   
%--------------------------------------------------------------------------

function eps_pl = evaluate_at_quad_points(hmsh, hspace, eps_pl_control_var)
    eps_pl = cell(hmsh.nlevels,1);
    
    ndofs_u = 0;
    for ilev = 1:hmsh.nlevels
        ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
        if (hmsh.nel_per_level(ilev) > 0)
            spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev});
            nsh_scalar = spu_lev.nsh_max/spu_lev.ncomp;
            sh_fun_u = reshape (spu_lev.shape_functions(1,:,1:nsh_scalar,:),  [hmsh.msh_lev{ilev}.nqn,  nsh_scalar, hmsh.msh_lev{ilev}.nel]);
    
            %dofs_u = 1:ndofs_u;
            %scalar_comp = 1:(hspace.ndof/hspace.ncomp);
    
            eps_pl_lev = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
            for iel=1:hmsh.msh_lev{ilev}.nel
                dofs_elem = spu_lev.connectivity(1:nsh_scalar,iel);
                eps_pl_lev(iel,:,:) = sh_fun_u(:,:,iel) * eps_pl_control_var(dofs_elem,:);
            end
    
            eps_pl{ilev} = eps_pl_lev;
    
        end
    end

end





%--------------------------------------------------------------------------
% subroutines   old
%--------------------------------------------------------------------------




function eps_pl = eps_pl_quasi_interpolant(hmsh_c, hmsh, hspace, eps_pl)

data = zeros(hmsh_c.mesh_of_level(1).nqn*hmsh_c.nel, 8); % 2D only!!!
id_lev =1;
for ilev = 1:hmsh_c.nlevels
    if (hmsh_c.nel_per_level(ilev) > 0)
    ngp_lev = hmsh_c.mesh_of_level(ilev).nqn * hmsh_c.nel_per_level(ilev);
    for idim = 1:hmsh_c.rdim        
        elements = hmsh_c.active{ilev};
        msh_lev = msh_evaluate_element_list( hmsh_c.mesh_of_level(ilev), elements );
        tmp =reshape (msh_lev.quad_nodes(idim,:,:), [ngp_lev,1]);
        data(id_lev:id_lev+ngp_lev-1,idim) = tmp;
    end

    data(id_lev:id_lev+ngp_lev-1,hmsh_c.rdim+1:end) = reshape(eps_pl{ilev}, [ngp_lev, 6]);
    id_lev = id_lev + ngp_lev;
    end
end


% Computing the coefficients of the quasi-interpolant in the space
% hspace based on the mesh hmsh, approximating the data
eps_pl_control_var = get_QI_coeffs(hspace,hmsh,data);

eps_pl = cell(hmsh.nlevels,1);

ndofs_u = 0;
for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
        spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev});
        nsh_scalar = spu_lev.nsh_max/spu_lev.ncomp;
        sh_fun_u = reshape (spu_lev.shape_functions(1,:,1:nsh_scalar,:),  [hmsh.msh_lev{ilev}.nqn,  nsh_scalar, hmsh.msh_lev{ilev}.nel]);

        %dofs_u = 1:ndofs_u;
        %scalar_comp = 1:(hspace.ndof/hspace.ncomp);

        eps_pl_lev = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
        for iel=1:hmsh.msh_lev{ilev}.nel
            dofs_elem = spu_lev.connectivity(1:nsh_scalar,iel);
            eps_pl_lev(iel,:,:) = sh_fun_u(:,:,iel) * eps_pl_control_var(dofs_elem,:);
        end

        eps_pl{ilev} = eps_pl_lev;

    end
end


end

function eps_pl = eps_pl_project_l2( hmsh_c, hspace_c, eps_pl_c, hmsh, hspace, Cref)

M = op_u_v_hier( hspace_c, hspace_c, hmsh_c);
rhs = eval_rhs_l2 (hspace_c, hmsh_c, eps_pl_c);

eps_pl_control_var=zeros(hspace.ndof/hspace.ncomp,6);
scalar_comp = 1:(hspace.ndof/hspace.ncomp);
scalar_comp_c = 1:(hspace_c.ndof/hspace_c.ncomp);
for i=1:6
    eps_pl_control_var_c = M(scalar_comp_c, scalar_comp_c)\rhs(:,i);
    eps_pl_control_var(:,i) = Cref(scalar_comp,scalar_comp_c) * eps_pl_control_var_c;
end

eps_pl = cell(hmsh.nlevels,1);

ndofs_u = 0;
for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
        spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev});
        nsh_scalar = spu_lev.nsh_max/spu_lev.ncomp;
        sh_fun_u = reshape (spu_lev.shape_functions(1,:,1:nsh_scalar,:),  [hmsh.msh_lev{ilev}.nqn,  nsh_scalar, hmsh.msh_lev{ilev}.nel]);

        %dofs_u = 1:ndofs_u;
        %scalar_comp = 1:(hspace.ndof/hspace.ncomp);

        eps_pl_lev = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
        for iel=1:hmsh.msh_lev{ilev}.nel
            dofs_elem = spu_lev.connectivity(1:nsh_scalar,iel);
            eps_pl_lev(iel,:,:) = sh_fun_u(:,:,iel) * eps_pl_control_var(dofs_elem,:);
        end

        eps_pl{ilev} = eps_pl_lev;

    end
end

end

function rhs = eval_rhs_l2 (hspv, hmsh, eps_pl)

rhs = zeros (hspv.ndof/hspv.ncomp, 6);
ndofs_v = 0;
for ilev = 1:hmsh.nlevels
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)

        spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), hmsh.msh_lev{ilev});
        rhs_lev = eval_rhs_lev (spv_lev, hmsh.msh_lev{ilev}, eps_pl{ilev}(:,:,:));

        dofs_v = 1:ndofs_v;
        scalar_comp = 1:(hspv.ndof/hspv.ncomp);
        rhs(dofs_v(scalar_comp),:) = rhs(dofs_v(scalar_comp),:) + hspv.Csub{ilev}(scalar_comp,scalar_comp)'*rhs_lev;
    end
end
end

function rhs = eval_rhs_lev (spv, msh, eps_pl)
nsh_scalar = spv.nsh_max/spv.ncomp;
v = reshape (spv.shape_functions(1,:,1:nsh_scalar,:),  [msh.nqn,  nsh_scalar, msh.nel]);
rhs = zeros(spv.ndof/spv.ncomp,6);
jacdet_weights = msh.jacdet .* msh.quad_weights;

for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))

        eps_pl_iel = reshape(eps_pl(iel,:,:), [msh.nqn,6])';
        rhs_iel = zeros(spv.nsh(iel)/spv.ncomp, 6);

        for igp =1: msh.nqn
            eps_pl_igp = eps_pl_iel(:,igp);
            rhs_iel = rhs_iel + (v(igp,:,iel).* jacdet_weights(igp,iel))' * eps_pl_igp';
        end

        % assembly vector
        rhs(spv.connectivity(1:nsh_scalar,iel),:)= rhs(spv.connectivity(1:nsh_scalar,iel),:) + rhs_iel;

    else
        warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
end

end