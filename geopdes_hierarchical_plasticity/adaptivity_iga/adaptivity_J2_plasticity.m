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

function [geometry, cell_hmsh, cell_hspace,  cell_hspace_scalar,  cell_u, cell_eps_pl, cell_sigma, solution_data] = adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data, plot_data)

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
cell_hspace_scalar = cell(method_data.nload,1);
cell_u = cell(method_data.nload,1);
cell_eps_pl = cell(method_data.nload,1); 
cell_sigma = cell(method_data.nload,1); 

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
    fprintf('----------------------------------------------------- Load step %d -----------------------------------------------------\n',iLoad);
    % ADAPTIVE LOOP
    iter = 0;
    while (1)
        iter = iter + 1;

        if (plot_data.print_info)
            fprintf('%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%\n',iter);
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

        sigma = cell(hmsh.nlevels,1);
        for ilev = 1:hmsh.nlevels
            sigma{ilev} = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
        end

        [u, eps_pl, sigma] = solve_J2_plasticity_hier (problem_data, method_data, adaptivity_data, hspace, hmsh, iLoad, u, eps_pl, sigma);



        % interpolate sigma
        sigma_store = zeros(hspace_scalar.ndof, 6); % control variables
        sigma_store(:,:) = history_variable_projection_hier(hspace_scalar, hmsh, sigma,  method_data.type_projection);

        % QI or L2 for eps_pl projection      
        eps_pl_store = zeros(hspace_scalar.ndof, 6); % control variables
        eps_pl_store(:,:) = history_variable_projection_hier(hspace_scalar, hmsh, eps_pl,  method_data.type_projection);



        nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;

       

        % ESTIMATE
        if (plot_data.print_info); disp('ESTIMATE:'); end
        %est = ones(hmsh.nel,1);% adaptivity_estimate_linear_el (u, hmsh, hspace, problem_data, adaptivity_data);
        est =  adaptivity_estimate_div_sigma_el (sigma_store, hmsh, hspace, hspace_scalar, problem_data, adaptivity_data);
        est_elem{iter} = est.';
        gest(iter) = norm (est);
        if (plot_data.print_info); fprintf('Computed error estimate: %e \n', gest(iter)); end

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
        

        if numel(marked) >= adaptivity_data.max_level -1            
            marked{adaptivity_data.max_level-1} = double.empty(0,1);
            num_marked =0;
            for ilev_mark =1:numel(marked)
                num_marked = num_marked + numel(marked{ilev_mark});
            end            
        end
        if num_marked == 0
            disp('Warning: no element refined')
        end
        
        if (plot_data.print_info)
            fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
            disp('REFINE:')
        end
        
        % REFINE
        hmsh_coarse = hmsh;
        hspace_coarse = hspace;
        hspace_scalar_coarse  = hspace_scalar;
        [hmsh, hspace, Cref] = adaptivity_refine (hmsh_coarse, hspace_coarse, marked, adaptivity_data);
        [~, hspace_scalar, Cref_scalar] = adaptivity_refine (hmsh_coarse, hspace_scalar_coarse, marked, adaptivity_data);

        % refine variables
        u = Cref * u;
        eps_pl_store = Cref_scalar*eps_pl_store; % control variables
        sigma_store = Cref_scalar*sigma_store; % control variables
        eps_pl = evaluate_at_quad_points(hmsh, hspace_scalar, eps_pl_store);



        fprintf('\n');
    end % end adaptivity step

    % coarsening step
    hmsh_fine = hmsh;
    [hmsh, hspace, u] =coarsening( hspace, hmsh_fine, u, est, adaptivity_data);
    [hmsh, hspace_scalar, tmp] =coarsening( hspace_scalar, hmsh_fine, [eps_pl_store, sigma_store] , est, adaptivity_data);
    eps_pl_store = tmp(:,1:size(eps_pl_store,2));
    eps_pl = evaluate_at_quad_points(hmsh, hspace_scalar, eps_pl_store);
    sigma_store = tmp(:,size(eps_pl_store,2)+1:end);


    % store data
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
    cell_hspace_scalar{iLoad} = hspace_scalar;
    cell_u{iLoad} = u;
    cell_eps_pl{iLoad} = eps_pl_store;
    cell_sigma{iLoad} = sigma_store;

    % plot hierarchical mesh
    % figure(1000)
    % hmsh_plot_cells (hmsh)
    % view(0,90)
    % 
    % drawnow


end % end load step

end





%--------------------------------------------------------------------------
% subroutines   
%--------------------------------------------------------------------------

function eps_pl = evaluate_at_quad_points(hmsh, hspace, eps_pl_control_var)
    eps_pl = cell(hmsh.nlevels,1);
    
    ndofs_u = 0;
    last_dof = cumsum (hspace.ndof_per_level);

    for ilev = 1:hmsh.nlevels
        ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
        if (hmsh.nel_per_level(ilev) > 0)
            spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev});
            % spu_lev = change_connectivity_localized_Csub (spu_lev, hspace, ilev);
            % nsh_scalar = spu_lev.nsh_max;
            sh_fun_u = spu_lev.shape_functions;
            %dofs_u = 1:ndofs_u;
            %scalar_comp = 1:(hspace.ndof/hspace.ncomp);    
            eps_pl_lev = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
            eps_pl_control_lev = hspace.Csub{ilev} * eps_pl_control_var(1:last_dof(ilev),:);
            for iel=1:hmsh.msh_lev{ilev}.nel               
                dofs_elem = spu_lev.connectivity(:,iel);
                eps_pl_lev(iel,:,:) = sh_fun_u(:,:,iel) * eps_pl_control_lev(dofs_elem,:);
            end
    
            eps_pl{ilev} = eps_pl_lev;
    
        end
    end

end


function [hmsh, hspace, u] =coarsening( hspace, hmsh, u, est, adaptivity_data)

 %% COARSENING =============================================================
        % MARK COARSENING
        disp('MARK COARSENING:')
        [marked_coarse, num_marked_coarse] = adaptivity_mark_coarsening (est, hmsh, hspace, adaptivity_data);

        % coarse only after the first time step if it also refines
        % COARSE
        if ~isempty(marked_coarse)
            
            fprintf('%d %s marked for coarsening \n', num_marked_coarse, adaptivity_data.flag);

            % Project the previous solution mesh onto the next refined mesh
            [hmsh_coarse, hspace_coarse, C_coar] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);
            u = C_coar * u;
            hmsh = hmsh_coarse;
            hspace = hspace_coarse;
        end
       
end

  







