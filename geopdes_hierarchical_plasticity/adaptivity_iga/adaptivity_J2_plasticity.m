% ADAPTIVITY_J2_PLASTICITY: adaptive hierarchical isogeometric J2 plasticity with projected history variables
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

function [geometry, cell_hmsh, cell_hspace,  cell_hspace_scalar,  cell_u, cell_eps_pl, cell_sigma, solution_data, cell_hmsh_scalar] = adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data, plot_data)

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

if (plot_data.plot_hmesh)
    fig_mesh = figure;
end
if (plot_data.plot_discrete_sol)
    fig_sol = figure;
end
nel = zeros (1, adaptivity_data.num_max_iter); ndof = nel; gest = nel+1;

cell_hmsh = cell(method_data.nload,1);
cell_hspace = cell(method_data.nload,1);
cell_hmsh_scalar = cell(method_data.nload,1);
cell_hspace_scalar = cell(method_data.nload,1);
cell_u = cell(method_data.nload,1);
cell_eps_pl = cell(method_data.nload,1); 
cell_sigma = cell(method_data.nload,1); 

[hmsh, hspace, geometry] = adaptivity_initialize_vector (problem_data, method_data);
hmsh_scalar = hmsh;


if isfield(method_data, 'degree_projection')
    degree_projection = method_data.degree_projection;
    regularity_projection = degree_projection - 1;        % max regularity (C^{p-1})
else
    degree_projection = method_data.degree;
    regularity_projection = method_data.regularity;
end
if any (strcmpi(method_data.type_projection, {'QI_C0', 'BEZIER_C0'}))
    regularity_projection = zeros(size(degree_projection)); % C^0
elseif strcmpi(method_data.type_projection, 'QI_graded')
end

[knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub_coarse-1, degree_projection, regularity_projection);
rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, geometry);
space_scalar =sp_bspline (knots, degree_projection, msh);
clear knots zeta rule qn  qw msh
hspace_scalar   = hierarchical_space (hmsh_scalar, space_scalar, method_data.space_type, method_data.truncated, regularity_projection);
hspace_dummy = hspace_scalar;

% QI_graded: C^{p-1} up to graded_transition_level, C^0 on finer levels
if strcmpi(method_data.type_projection, 'QI_graded')
    if isfield(method_data, 'graded_transition_level')
        graded_trans = method_data.graded_transition_level;
    else
        graded_trans = 2;   % default: levels 1-2 smooth, 3+ are C^0
    end
    max_reg = degree_projection - 1;
    min_reg = zeros(size(degree_projection));
    n_alloc = max(adaptivity_data.max_level + 2, 10);
    reg_per_level = cell(1, n_alloc);
    for lev = 1:n_alloc
        if lev <= graded_trans
            reg_per_level{lev} = max_reg;
        else
            reg_per_level{lev} = min_reg;
        end
    end
    hspace_scalar.regularity_per_level = reg_per_level;

    ndim = numel(degree_projection);
    if ~isempty(hspace_scalar.boundary)
        for iside = 1:numel(hspace_scalar.boundary)
            ind = setdiff(1:ndim, ceil(iside/2));
            bnd_reg_per_level = cell(1, n_alloc);
            for lev = 1:n_alloc
                if lev <= graded_trans
                    bnd_reg_per_level{lev} = max_reg(ind);
                else
                    bnd_reg_per_level{lev} = min_reg(ind);
                end
            end
            hspace_scalar.boundary(iside).regularity_per_level = bnd_reg_per_level;
        end
    end

    fprintf('  QI_graded: levels 1-%d use C^%d, levels %d+ use C^0\n', ...
        graded_trans, max_reg(1), graded_trans+1);
end

if strcmpi(method_data.type_projection, 'QI_graded')
    projection_type = 'QI';
else
    projection_type = method_data.type_projection;
end


if isfield(method_data, 'num_bisections')
    num_bisections = method_data.num_bisections;
else
    num_bisections = 0;
end

if isfield(method_data, 'qi_lambda')
    qi_lambda = method_data.qi_lambda;
else
    qi_lambda = [];
end
if any (strcmpi (method_data.type_projection, {'BEZIER', 'BEZIER_C0', 'DLSQ', 'DLSQ_W'})) && num_bisections > 0
    error ('adaptivity_J2_plasticity: type_projection ''%s'' is incompatible with num_bisections > 0', method_data.type_projection);
end

% CPT and WPLSQ transfer quadrature data directly (Hennig et al. 2018), PROJECT round-trips through the coefficients
if isfield(method_data, 'type_transfer')
    type_transfer = method_data.type_transfer;
else
    type_transfer = 'PROJECT';
end
if ~any (strcmpi (type_transfer, {'PROJECT', 'PROJECT_NEW', 'CPT', 'WPLSQ'}))
    error ('adaptivity_J2_plasticity: unknown type_transfer ''%s''', type_transfer);
end
if any (strcmpi (type_transfer, {'PROJECT_NEW', 'CPT', 'WPLSQ'})) && num_bisections > 0
    error ('adaptivity_J2_plasticity: type_transfer ''%s'' is incompatible with num_bisections > 0', type_transfer);
end
if num_bisections > 0
    [hmsh_scalar, hspace_scalar] = refine_projection_space(hmsh_scalar, hspace_scalar, adaptivity_data, num_bisections);
end


u = zeros (hspace.ndof,1);
eps_pl = cell(hmsh.nlevels,1);
for ilev = 1:hmsh.nlevels
    eps_pl{ilev} = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
end


if (isfield (adaptivity_data, 'max_level') && ...
    adaptivity_data.num_max_iter < adaptivity_data.max_level)
  warning ('adaptivity_J2_plasticity:tooFewIterations', ...
           ['num_max_iter = %d gives %d refinement passes per load step, which ', ...
            'cannot lift an element from level 1 to max_level = %d. The mesh may ', ...
            'lag behind a moving front. Use num_max_iter >= max_level.'], ...
           adaptivity_data.num_max_iter, adaptivity_data.num_max_iter - 1, ...
           adaptivity_data.max_level);
end

for iLoad = 1: method_data.nload
    fprintf('----------------------------------------------------- Load step %d -----------------------------------------------------\n',iLoad);
    iter = 0;
    est = [];   % may stay empty when num_max_iter == 1 (no estimator computed)
    while (1)
        iter = iter + 1;

        if (plot_data.print_info)
            fprintf('%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%\n',iter);
        end
        if (~hspace_check_partition_of_unity (hspace, hmsh))
            disp('ERROR: The partition-of-the-unity property does not hold.')
            solution_data.flag = -1; break
        end

        if (plot_data.print_info)
            disp('SOLVE:')
            fprintf('Number of elements: %d. Total DOFs: %d \n', hmsh.nel, hspace.ndof);
        end

        sigma = cell(hmsh.nlevels,1);
        for ilev = 1:hmsh.nlevels
            sigma{ilev} = zeros(hmsh.nel_per_level(ilev),hmsh.mesh_of_level(ilev).nqn,6);
        end

        [u, eps_pl, sigma] = solve_J2_plasticity_hier (problem_data, method_data, adaptivity_data, hspace, hmsh, iLoad, u, eps_pl, sigma);

        if (plot_data.print_info); disp('PROJECT STATE VARIABLES:'); end
        % project sigma and eps_pl as one 12-component field to share the per-level precompute
        hist_both = cell(size(sigma));
        for ilev = 1:numel(sigma)
            if ~isempty(sigma{ilev})
                hist_both{ilev} = cat(3, sigma{ilev}, eps_pl{ilev});
            else
                hist_both{ilev} = sigma{ilev};
            end
        end
        both_store   = history_variable_projection_hier(hspace_scalar, hmsh_scalar, hist_both, projection_type, hmsh, qi_lambda);
        sigma_store  = both_store(:, 1:6);
        eps_pl_store = both_store(:, 7:12);

        % Opt-in round-trip diagnostic, a ratio above one means the projection inflated the field
        global PROJ_DIAG %#ok<GVMIS>
        if (~isempty(PROJ_DIAG) && PROJ_DIAG)
            back = evaluate_at_quad_points(hmsh_scalar, hspace_scalar, sigma_store);
            din = 0; dout = 0; dmax = 0;
            for il = 1:numel(sigma)
                if (isempty(sigma{il})); continue; end
                din  = max(din,  max(abs(sigma{il}(:))));
                dout = max(dout, max(abs(back{il}(:))));
                d    = abs(back{il}(:) - sigma{il}(:));
                dmax = max(dmax, max(d));
            end
            fprintf(['[PROJ_DIAG] load %d iter %d  nel %d  max|in| %.4e  ' ...
                     'max|out| %.4e  ratio %.4f  max|out-in| %.4e  max|coef| %.4e\n'], ...
                    iLoad, iter, hmsh.nel, din, dout, dout/max(din,realmin), dmax, ...
                    max(abs(both_store(:))));
        end

        nel(iter) = hmsh.nel; ndof(iter) = hspace.ndof;
        
        if (iter > adaptivity_data.num_max_iter)
            disp('skip refinement')
            solution_data.flag = 2; break
        end


        if (plot_data.print_info); disp('ESTIMATE:'); end
        if (~isfield(adaptivity_data, 'estimator'))
            adaptivity_data.estimator = 'div_sigma_direct';
        end

        switch lower(adaptivity_data.estimator)
            case {'div_sigma_direct', 'div_sigma_direct_bdry'}
                % The direct estimator pairs the scalar space level-by-level with the primal
                % mesh, so remap onto hspace_dummy when the projection mesh is bisected
                if num_bisections > 0
                    eps_pl_store_est = zeros(hspace_dummy.ndof, 6);
                    eps_pl_store_est(:,:) = history_variable_projection_hier(hspace_dummy, hmsh, eps_pl, 'L2', hmsh);
                    sigma_store_est = zeros(hspace_dummy.ndof, 6);
                    sigma_store_est(:,:) = history_variable_projection_hier(hspace_dummy, hmsh, sigma, 'L2', hmsh);
                    hspace_est = hspace_dummy;
                else
                    eps_pl_store_est = eps_pl_store;
                    sigma_store_est  = sigma_store;
                    hspace_est = hspace_scalar;
                end

                if strcmpi(adaptivity_data.estimator, 'div_sigma_direct_bdry')
                    adaptivity_data.boundary_term_full = true;
                    adaptivity_data.load_mult = iLoad / method_data.nload;
                end

                est = adaptivity_estimate_div_sigma_direct_el (u, eps_pl_store_est, ...
                    sigma_store_est, geometry, hmsh, hspace, hspace_est, ...
                    problem_data, adaptivity_data);

            case {'plastic_front'}
                % Solution-independent marker, every projection variant gets the same mesh sequence
                est = adaptivity_estimate_distance_plastic_front (hmsh, ...
                    problem_data, iLoad / method_data.nload);

            case {'plastic_front_sphere'}
                est = adaptivity_estimate_distance_plastic_front_sphere (hmsh, ...
                    problem_data, iLoad / method_data.nload);

            case {'plastic_front_discrete'}
                % Front taken as the boundary of {sigma_vm >= sigma_y} in the computed stress
                if num_bisections > 0
                    sigma_store_est = zeros(hspace_dummy.ndof, 6);
                    sigma_store_est(:,:) = history_variable_projection_hier(hspace_dummy, hmsh, sigma, 'L2', hmsh);
                    hspace_est = hspace_dummy;
                else
                    sigma_store_est = sigma_store;
                    hspace_est = hspace_scalar;
                end
                est = adaptivity_estimate_plastic_front_discrete_el (...
                    sigma_store_est, hmsh, hspace_est, problem_data, adaptivity_data);

            case {'div_sigma', 'residual', 'stress_gradient', 'jump', ...
                   'plastic_front_spere', 'sphere_front'}
                if num_bisections > 0
                    sigma_store_est = zeros(hspace_dummy.ndof, 6);
                    sigma_store_est(:,:) = history_variable_projection_hier(hspace_dummy, hmsh, sigma, 'L2', hmsh);
                    hspace_est = hspace_dummy;
                    hmsh_est   = hmsh;
                else
                    sigma_store_est = sigma_store;
                    hspace_est = hspace_scalar;
                    hmsh_est   = hmsh;
                end

                switch lower(adaptivity_data.estimator)
                    case {'div_sigma', 'residual'}
                        est = adaptivity_estimate_div_sigma_el (sigma_store_est, geometry, hmsh, hspace, hmsh_est, hspace_est, problem_data, adaptivity_data);
                    case {'stress_gradient', 'jump'}
                        est = adaptivity_estimate_stress_gradient_el (sigma_store_est, geometry, hmsh, hspace, hmsh_est, hspace_est, problem_data, adaptivity_data);
                    case {'plastic_front_spere', 'sphere_front'}
                        est = adaptivity_estimate_stress_gradient_el (sigma_store_est, geometry, hmsh, hspace, hmsh_est, hspace_est, problem_data, adaptivity_data);
                end

            otherwise
                error('adaptivity_J2_plasticity:estimator', ...
                      'Unknown adaptivity_data.estimator = %s', adaptivity_data.estimator);
        end
        est_elem{iter} = est.';
        gest(iter) = norm (est);
        if (plot_data.print_info); fprintf('Computed error estimate: %e \n', gest(iter)); end

        if (gest(iter) < adaptivity_data.tol)
            disp('Success: The error estimation reached the desired tolerance');
            solution_data.flag = 1; break
        elseif (iter == adaptivity_data.num_max_iter)
            disp('Warning: reached the maximum number of iterations')
            solution_data.flag = 2; break
        % Deliberately no global nlevels >= max_level stop here, the per-level mark cap
        % below enforces max_level element-wise without freezing the whole mesh
        elseif (hspace.ndof > adaptivity_data.max_ndof)
            disp('Warning: reached the maximum number of DOFs')
            solution_data.flag = 4; break
        elseif (hmsh.nel > adaptivity_data.max_nel)
            disp('Warning: reached the maximum number of elements')
            solution_data.flag = 5; break
        end

        if (plot_data.print_info); disp('MARK:'); end
        [marked, num_marked] = adaptivity_mark (est, hmsh, hspace, adaptivity_data);


        % Drop marks at levels >= max_level so refinement never creates cells above it
        if numel(marked) >= adaptivity_data.max_level
            for k = adaptivity_data.max_level:numel(marked)
                marked{k} = double.empty(0,1);
            end
            num_marked = 0;
            for ilev_mark = 1:numel(marked)
                num_marked = num_marked + numel(marked{ilev_mark});
            end
        end
        if num_marked == 0
            disp('Warning: no element refined')
            break
        end

        % Admissibility expansion (Buffa and Giannelli, M3AS 2016) computed once on the
        % displacement hierarchy and reused by every refine call so all spaces share one mesh sequence
        num_marked_in = num_marked;
        marked = mark_admissible (hmsh, hspace, marked, adaptivity_data);
        num_marked = sum (cellfun (@numel, marked));

        if (plot_data.print_info)
            if (num_marked > num_marked_in)
                fprintf('%d %s marked for refinement (%d after admissibility class %d) \n', ...
                        num_marked_in, adaptivity_data.flag, num_marked, adaptivity_data.adm);
            else
                fprintf('%d %s marked for refinement \n', num_marked, adaptivity_data.flag);
            end
            disp('REFINE:')
        end

        hmsh_coarse = hmsh;
        hmsh_scalar_coarse = hmsh_scalar;
        hspace_coarse = hspace;
        hspace_scalar_coarse  = hspace_scalar;
        [hmsh, hspace, Cref] = adaptivity_refine (hmsh_coarse, hspace_coarse, marked, adaptivity_data);
        
        
        if num_bisections > 0
            [~, hspace_dummy] = adaptivity_refine (hmsh_coarse, hspace_dummy, marked, adaptivity_data);
             [hmsh_scalar, hspace_scalar] = refine_projection_space(hmsh, hspace_dummy, adaptivity_data, num_bisections);

             tmp_M = op_u_v_hier(hspace_scalar,hspace_scalar,hmsh_scalar);
             hspace_scalar_in_finer_mesh = hspace_in_finer_mesh(hspace_scalar_coarse, hmsh_scalar_coarse, hmsh_scalar);
             tmp_G = op_u_v_hier(hspace_scalar_in_finer_mesh,hspace_scalar,hmsh_scalar);
             tmp_lhs_e = zeros(hspace_scalar.ndof,6);
             tmp_lhs_s = zeros(hspace_scalar.ndof,6);
             for ivar = 1:6
                tmp_lhs_e(:,ivar) = tmp_M\ (tmp_G*eps_pl_store(:,ivar));
                tmp_lhs_s(:,ivar) = tmp_M\(tmp_G*sigma_store(:,ivar));
             end
             eps_pl_store = tmp_lhs_e;
             sigma_store =tmp_lhs_s;

        else
            [hmsh_scalar, hspace_scalar, Cref_scalar] = adaptivity_refine (hmsh_scalar_coarse, hspace_scalar_coarse, marked, adaptivity_data);
            eps_pl_store = Cref_scalar*eps_pl_store; % control variables
            sigma_store = Cref_scalar*sigma_store; % control variables


        end


        u = Cref * u;
        switch upper(type_transfer)
            case 'PROJECT'
                eps_pl = evaluate_at_quad_points(hmsh, hspace_scalar, eps_pl_store);
            case 'PROJECT_NEW'
                % Hennig et al. 2018, Eq. (40): write only the newly activated elements
                eps_pl = keep_unchanged_elements(hmsh_coarse, eps_pl, hmsh, evaluate_at_quad_points(hmsh, hspace_scalar, eps_pl_store));
            case 'CPT'
                eps_pl = cpt_transfer_hier(hmsh_coarse, eps_pl, hmsh);
            case 'WPLSQ'
                eps_pl = wplsq_transfer_hier(hmsh_coarse, hspace_scalar_coarse, eps_pl, hmsh);
        end



        fprintf('\n');
        
    end % end adaptivity step
    if (adaptivity_data.mark_param_coarsening==0. || isempty(est))
            disp('skip coarsening')
    else
        hmsh_fine = hmsh;
        % pre-coarsen scalar space, CPT and WPLSQ need it after hspace_scalar is overwritten
        hspace_scalar_fine_pre = hspace_scalar;
        [hmsh, hspace, u, reactivated_elements] =coarsening( hspace, hmsh_fine, u, est, adaptivity_data);

        if (hmsh.nel == hmsh_fine.nel)
        % nothing reactivated, skip the lossy projection round-trip
        elseif num_bisections > 0
        % Reuse the reactivated elements so hspace_dummy cannot drift away from hmsh
        [~, hspace_dummy] =coarsening( hspace_dummy, hmsh_fine, zeros(hspace_dummy.ndof,1), est, adaptivity_data, reactivated_elements);

        hmsh_scalar_fine = hmsh_scalar;
        hspace_scalar_fine = hspace_scalar;
        [hmsh_scalar, hspace_scalar] = refine_projection_space(hmsh, hspace_dummy, adaptivity_data, num_bisections);

        tmp_M = op_u_v_hier(hspace_scalar,hspace_scalar,hmsh_scalar);
        hspace_scalar_in_finer_mesh = hspace_in_finer_mesh(hspace_scalar, hmsh_scalar, hmsh_scalar_fine);

        tmp_G = op_u_v_hier(hspace_scalar_fine,hspace_scalar_in_finer_mesh,hmsh_scalar_fine);
        tmp_lhs_e = zeros(hspace_scalar.ndof,6);
        tmp_lhs_s = zeros(hspace_scalar.ndof,6);
        for ivar = 1:6
            tmp_lhs_e(:,ivar) = tmp_M\ (tmp_G*eps_pl_store(:,ivar));
            tmp_lhs_s(:,ivar) = tmp_M\(tmp_G*sigma_store(:,ivar));
        end
        eps_pl_store = tmp_lhs_e;
        sigma_store =tmp_lhs_s;

        else
            [hmsh_scalar, hspace_scalar, tmp] =coarsening( hspace_scalar, hmsh_fine, [eps_pl_store, sigma_store] , est, adaptivity_data, reactivated_elements);
            eps_pl_store = tmp(:,1:size(eps_pl_store,2));
            sigma_store = tmp(:,size(eps_pl_store,2)+1:end);

        end
        if (hmsh.nel ~= hmsh_fine.nel)
            switch upper(type_transfer)
                case 'PROJECT'
                    eps_pl = evaluate_at_quad_points(hmsh, hspace_scalar, eps_pl_store);
                case 'PROJECT_NEW'
                    % on coarsening the projection supplies all coarsened cells, survivors keep their data
                    eps_pl = keep_unchanged_elements(hmsh_fine, eps_pl, hmsh, evaluate_at_quad_points(hmsh, hspace_scalar, eps_pl_store));
                case 'CPT'
                    eps_pl = cpt_transfer_hier(hmsh_fine, eps_pl, hmsh);
                case 'WPLSQ'
                    eps_pl = wplsq_transfer_hier(hmsh_fine, hspace_scalar_fine_pre, eps_pl, hmsh);
            end
        end
    end

    


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


    cell_hmsh{iLoad} = hmsh;
    cell_hspace{iLoad} = hspace;
    cell_hspace_scalar{iLoad} = hspace_scalar;
    cell_hmsh_scalar{iLoad} = hmsh_scalar;
    cell_u{iLoad} = u;
    cell_eps_pl{iLoad} = eps_pl_store;
    cell_sigma{iLoad} = sigma_store;



end % end load step

end






function marked = mark_admissible (hmsh, hspace, marked, adaptivity_data)
% MARK_ADMISSIBLE: expand marks so refinement yields a mesh admissible of class adaptivity_data.adm
% hrefine returns a superset of the input marks, so the result feeds adaptivity_refine directly

    if (~isfield (adaptivity_data, 'adm_strategy') || isempty (adaptivity_data.adm_strategy))
        return
    end

    switch (lower (adaptivity_data.adm_strategy))
        case {'none', 'off'}
            % refinement is left unconstrained on purpose

        case 'admissible'
            if (~isfield (adaptivity_data, 'adm') || isempty (adaptivity_data.adm))
                error ('adaptivity_J2_plasticity:adm_missing', ...
                       ['adm_strategy = ''admissible'' requires the admissibility ' ...
                        'class in adaptivity_data.adm']);
            end
            % hrefine errors for m < 2, class 1 imposes no constraint anyway
            if (adaptivity_data.adm > 1)
                shape_in = size (marked);
                marked = hrefine (hmsh, hspace, marked, adaptivity_data.adm);
                % hrefine returns a row cell, restore the caller's orientation
                if (shape_in(1) >= shape_in(2))
                    marked = marked(:);
                else
                    marked = marked(:).';
                end
            end

        otherwise
            error ('adaptivity_J2_plasticity:adm_strategy_unknown', ...
                   ['unknown adm_strategy ''%s'' (supported: ''admissible'', ' ...
                    '''none''). Note adaptivity_refine_fsb, which the poisson ' ...
                    'and thermomech solvers dispatch to for ''balancing'', is ' ...
                    'not present in this tree.'], adaptivity_data.adm_strategy);
    end
end


function  [hmsh_scalar, hspace_scalar] = refine_projection_space(hmsh_scalar, hspace_scalar, adaptivity_data, num_bisections)
    for i=1: num_bisections
        marked = mark_bisect_mesh(hmsh_scalar);
        [hmsh_scalar, hspace_scalar] = adaptivity_refine (hmsh_scalar, hspace_scalar, marked, adaptivity_data);
    end
end


function rhs = op_Gu_hier (hspace, hmsh_fine, hspace_fine, uhat_fine)

  rhs = zeros (hspace.ndof, 1);

  ndofs = 0;
  for ilev = 1:hmsh_fine.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh_fine.nel_per_level(ilev) > 0)
      x = cell (hmsh_fine.rdim, 1);
      for idim = 1:hmsh_fine.rdim
        x{idim} = reshape (hmsh_fine.msh_lev{ilev}.geo_map(idim,:,:), hmsh_fine.msh_lev{ilev}.nqn, hmsh_fine.nel_per_level(ilev));
      end

      sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh_fine.msh_lev{ilev}, 'value', true);      
      sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);

      ndof_until_lev = sum (hspace_fine.ndof_per_level(1:ilev));
      uhat_lev = hspace_fine.Csub{ilev} * uhat_fine(1:ndof_until_lev);
      sp_lev_fine = sp_evaluate_element_list (hspace_fine.space_of_level(ilev), hmsh_fine.msh_lev{ilev}, 'value', true); 
      sp_lev_fine = change_connectivity_localized_Csub (sp_lev_fine, hspace_fine, ilev);
      utemp = sp_eval_msh (uhat_lev, sp_lev_fine, hmsh_fine.msh_lev{ilev}, {'value'});
      u_fine = utemp{1};

      b_lev = op_f_v (sp_lev, hmsh_fine.msh_lev{ilev}, u_fine);

      dofs = 1:ndofs;
      rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    end
  end

end



function eps_pl = evaluate_at_quad_points(hmsh, hspace, eps_pl_control_var)
    % Level-by-level pairing is only valid when hspace and hmsh share a hierarchy,
    % otherwise it silently returns zeros, so use the nested evaluator
    if (hspace.nlevels ~= hmsh.nlevels)
        eps_pl = hspace_eval_hmsh_nested (eps_pl_control_var, hspace, hmsh);
        return
    end

    eps_pl = cell(hmsh.nlevels,1);

    ndofs_u = 0;
    last_dof = cumsum (hspace.ndof_per_level);

    for ilev = 1:hmsh.nlevels
        ndofs_u = ndofs_u + hspace.ndof_per_level(ilev);
        if (hmsh.nel_per_level(ilev) > 0)
            spu_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev});
            sh_fun_u = spu_lev.shape_functions;
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


function [hmsh, hspace, u, reactivated_elements] =coarsening( hspace, hmsh, u, est, adaptivity_data, reactivated_elements_in)

        disp('MARK COARSENING:')
        [marked_coarse, num_marked_coarse] = adaptivity_mark_coarsening (est, hmsh, hspace, adaptivity_data);

        reactivated_elements = {};

        if ~isempty(marked_coarse)

            fprintf('%d %s marked for coarsening \n', num_marked_coarse, adaptivity_data.flag);

            if nargin >= 6
                [hmsh_coarse, hspace_coarse, C_coar, reactivated_elements] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data, reactivated_elements_in);
            else
                [hmsh_coarse, hspace_coarse, C_coar, reactivated_elements] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);
            end
            u = C_coar * u;
            hmsh = hmsh_coarse;
            hspace = hspace_coarse;
        end

end

  







% Hennig et al., CMAME 334 (2018), Eqs. (39) and (40): surviving elements keep their old
% quadrature data, only new elements take the projected values. A surviving cell keeps its
% (level, cell index) and quadrature points, so matching on them is exact
function eps_new = keep_unchanged_elements (hmsh_old, eps_old, hmsh_new, eps_proj)
  eps_new = eps_proj;
  for ilev = 1:min (hmsh_old.nlevels, hmsh_new.nlevels)
    if (numel (eps_old) < ilev); continue; end
    new_idx = hmsh_new.active{ilev}(:);
    old_idx = hmsh_old.active{ilev}(:);
    [kept, where] = ismember (new_idx, old_idx);
    if (~any (kept)); continue; end
    % guard against a stale eps_old that does not match its own mesh
    if (size (eps_old{ilev}, 1) ~= numel (old_idx) || size (eps_new{ilev}, 1) ~= numel (new_idx) ...
        || size (eps_old{ilev}, 2) ~= size (eps_new{ilev}, 2)); continue; end
    eps_new{ilev}(kept, :, :) = eps_old{ilev}(where(kept), :, :);
  end
end
