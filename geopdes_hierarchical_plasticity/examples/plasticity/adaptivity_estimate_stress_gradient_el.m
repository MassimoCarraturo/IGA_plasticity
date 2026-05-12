function est = adaptivity_estimate_stress_gradient_el(sigma_store, geometry, hmsh, hspace, hmsh_scalar, hspace_scalar, problem_data, adaptivity_data)
% ADAPTIVITY_ESTIMATE_STRESS_GRADIENT_EL  Element-local error indicator
% based on the stress-gradient seminorm (proxy for the inter-element
% stress jump):
%
%       eta_K = h_K * sqrt( integral_K  ||grad sigma||_F^2  dK )
%
% where ||grad sigma||_F is the Frobenius norm of the gradient tensor of
% the projected stress field (summed over all 6 Voigt components).
%
% Rationale: for a smooth basis the projected stress is C^0 (or higher),
% so true inter-element jumps are zero.  The stress *gradient* is the
% natural proxy: where ||grad sigma|| is large the field is changing
% rapidly within the element, and any piecewise-constant reconstruction
% would produce large face jumps there.  In the eighth-sphere benchmark
% the elastic-plastic front is exactly such a region; this indicator
% concentrates refinement there and leaves the (nearly affine) elastic
% interior coarse, matching the user's intent.
%
% Same input/output convention as adaptivity_estimate_div_sigma_el:
%   est is a vector of size [hmsh.nel x 1] when adaptivity_data.flag is
%   'elements', or [hspace.ndof x 1] when 'functions'.
%
% This estimator uses the bulk-evaluation routine hspace_eval_hmsh, which
% is several orders of magnitude faster than the point-by-point
% sp_eval_phys loop used by the original residual-divergence estimator.
%
% Caveat: hspace_eval_hmsh requires hspace_scalar to live on the same
% hierarchical mesh as hmsh.  This holds for the L2 / QI / QI_C0
% (num_bisections == 0) configurations used in the sphere benchmark.
% For a refined projection mesh (QI_C0 with num_bisections > 0), fall
% back to the slow per-point evaluation by setting
% adaptivity_data.gradient_estimator_legacy = true.

    if (isfield(adaptivity_data, 'C0_est'))
        C0_est = adaptivity_data.C0_est;
    else
        C0_est = 1;
    end
    use_legacy = isfield(adaptivity_data, 'gradient_estimator_legacy') ...
                 && adaptivity_data.gradient_estimator_legacy;

    ncomp_sigma = size(sigma_store, 2);
    ndim        = hmsh.ndim;
    nqn_total   = hmsh.mesh_of_level(1).nqn;

    % grad_sq(iquad, iel) = sum over Voigt components and physical-grad
    % directions of |dsigma_c / dx_d|^2
    grad_sq = zeros(nqn_total, hmsh.nel);

    if (~use_legacy)
        % Fast path: bulk hspace_eval_hmsh on the projection space.
        % Returns a [ndim x nqn x nel] gradient array per Voigt component.
        for ifield = 1:ncomp_sigma
            [hg, F] = hspace_eval_hmsh (sigma_store(:, ifield), ...
                                        hspace_scalar, hmsh_scalar, 'gradient');
            % hg has shape [ndim x nqn x nel]
            grad_sq = grad_sq + reshape (sum (hg.^2, 1), nqn_total, hmsh.nel);
        end
    else
        % Legacy slow path (e.g. for hmsh_scalar != hmsh).
        [~, F] = hspace_eval_hmsh (zeros(hspace.ndof, 1), hspace, hmsh, 'value');
        for ifield = 1:ncomp_sigma
            for iel = 1:hmsh.nel
                for iquad = 1:nqn_total
                    g = sp_eval_phys (sigma_store(:, ifield), hspace_scalar, ...
                                      hmsh_scalar, geometry, F(:, iquad, iel), ...
                                      'gradient');
                    grad_sq(iquad, iel) = grad_sq(iquad, iel) + sum(g(:).^2);
                end
            end
        end
    end

    switch adaptivity_data.flag
        case 'elements'
            w = []; h = [];
            for ilev = 1:hmsh.nlevels
                if (hmsh.msh_lev{ilev}.nel ~= 0)
                    w = cat (2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
                    h = cat (1, h, hmsh.msh_lev{ilev}.element_size(:));
                end
            end
            h = h * sqrt (hmsh.ndim);

            % Integrate ||grad sigma||_F^2 over each element
            elem_int = sum (grad_sq .* w, 1);

            est = (h(:)).^2 .* elem_int(:);
            est = C0_est * sqrt (est);

        case 'functions'
            % Function-based variant: distribute the elementwise indicator
            % to the supporting basis functions weighted by partition-of-unity
            ms = zeros (hmsh.nlevels, 1);
            for ilev = 1:hmsh.nlevels
                if (hmsh.msh_lev{ilev}.nel ~= 0)
                    ms(ilev) = max (hmsh.msh_lev{ilev}.element_size);
                else
                    ms(ilev) = 0;
                end
            end
            ms = ms * sqrt (hmsh.ndim);

            Nf = cumsum ([0; hspace.ndof_per_level(:)]);
            dof_level = zeros (hspace.ndof, 1);
            for lev = 1:hspace.nlevels
                dof_level(Nf(lev)+1:Nf(lev+1)) = lev;
            end
            coef = ms(dof_level).^2 .* hspace.coeff_pou(:);

            est = zeros (hspace.ndof, 1);
            ndofs = 0;
            Ne = cumsum ([0; hmsh.nel_per_level(:)]);
            for ilev = 1:hmsh.nlevels
                ndofs = ndofs + hspace.ndof_per_level(ilev);
                if (hmsh.nel_per_level(ilev) > 0)
                    ind_e = (Ne(ilev)+1):Ne(ilev+1);
                    sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
                    b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, ...
                                    reshape (grad_sq(:, ind_e), [1, nqn_total, numel(ind_e)]));
                    dofs = 1:ndofs;
                    est(dofs) = est(dofs) + hspace.Csub{ilev}.' * b_lev;
                end
            end
            est = coef .* est;
            est = C0_est * sqrt (est);

        otherwise
            error('adaptivity_estimate_stress_gradient_el:flag', ...
                  'Unknown adaptivity_data.flag = %s', adaptivity_data.flag);
    end
end
