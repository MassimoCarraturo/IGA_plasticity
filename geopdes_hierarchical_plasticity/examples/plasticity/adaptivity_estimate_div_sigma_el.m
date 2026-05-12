function est = adaptivity_estimate_div_sigma_el (sigma_store, geometry, hmsh, hspace, hmsh_scalar, hspace_scalar, problem_data, adaptivity_data)
% ADAPTIVITY_ESTIMATE_DIV_SIGMA_EL  Element-local residual error indicator
% based on the strong-form equilibrium residual:
%
%       eta_K = h_K * || f + div(sigma_h) ||_{L^2(K)}
%
% The residual f + div(sigma) is zero for the exact solution; its magnitude
% measures how well the discrete stress field locally satisfies equilibrium.
%
% Voigt convention: 1=sig_xx, 2=sig_yy, 3=sig_zz, 4=sig_xy, 5=sig_yz, 6=sig_xz
%
% Two evaluation paths:
%   Fast (default): bulk hspace_eval_hmsh on the projection space.
%     Requires hmsh_scalar to coincide with hmsh (true when
%     num_bisections == 0, i.e. L2, QI, QI_C0 with the default setup).
%   Legacy: per-point sp_eval_phys via Newton inverse-map.
%     Activated by setting adaptivity_data.div_sigma_estimator_legacy = true.
%
% Same input/output convention as adaptivity_estimate_stress_gradient_el.

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end
use_legacy = isfield(adaptivity_data, 'div_sigma_estimator_legacy') ...
             && adaptivity_data.div_sigma_estimator_legacy;

ncomp      = hspace.ncomp;          % 2 (2D) or 3 (3D)
ncomp_sigma = size(sigma_store, 2); % 6 Voigt components (or 3 in 2D: xx, yy, xy)
ndim       = hmsh.ndim;
nqn_total  = hmsh.mesh_of_level(1).nqn;

% Voigt-to-divergence map:  divergence(i) += d sigma_{ij} / d x_j
% For each Voigt index c, store which row of div it contributes to and
% which derivative direction(s) to take.
%   c  | sigma |  div(sigma)_row  | derivative dir
%   1  | xx    |       1          |     1  (d/dx)
%   2  | yy    |       2          |     2  (d/dy)
%   3  | zz    |       3          |     3  (d/dz)  [3D only]
%   4  | xy    |     1 and 2      |   2 and 1
%   5  | yz    |     2 and 3      |   3 and 2      [3D only]
%   6  | xz    |     1 and 3      |   3 and 1      [3D only]

divergence = zeros(ncomp, nqn_total, hmsh.nel);

if (~use_legacy)
    % ---- Fast path: bulk evaluation via hspace_eval_hmsh ----------------
    for ifield = 1:ncomp_sigma
        [hgrad, F] = hspace_eval_hmsh (sigma_store(:, ifield), ...
                                        hspace_scalar, hmsh_scalar, 'gradient');
        % hgrad has shape [ndim x nqn x nel]
        divergence = accumulate_divergence(divergence, hgrad, ifield, ncomp);
    end
else
    % ---- Legacy path: per-point sp_eval_phys ----------------------------
    [~, F] = hspace_eval_hmsh (zeros(hspace.ndof, 1), hspace, hmsh, 'value');
    for ifield = 1:ncomp_sigma
        hgrad = zeros(ndim, nqn_total, hmsh.nel);
        for iel = 1:hmsh.nel
            for iquad = 1:nqn_total
                tmp = sp_eval_phys (sigma_store(:, ifield), hspace_scalar, ...
                                    hmsh_scalar, geometry, F(:, iquad, iel), ...
                                    'gradient');
                hgrad(:, iquad, iel) = tmp;
            end
        end
        divergence = accumulate_divergence(divergence, hgrad, ifield, ncomp);
    end
end

% Evaluate the body force f at the physical quadrature points
x = cell (hmsh.rdim, 1);
for idim = 1:hmsh.rdim
    x{idim} = reshape (F(idim, :, :), [], hmsh.nel);
end
valf = problem_data.f (x{:});

% Residual: (f + div sigma)^2, summed over components
aux = sum((valf + divergence).^2, 1);
aux = reshape (aux, nqn_total, hmsh.nel);

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

        % Integrate ||f + div sigma||^2 over each element
        elem_int = sum (aux .* w, 1);
        est = (h(:)).^2 .* elem_int(:);
        est = C0_est * sqrt (est);

    case 'functions'
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

        est = zeros(hspace.ndof, 1);
        ndofs = 0;
        Ne = cumsum([0; hmsh.nel_per_level(:)]);
        for ilev = 1:hmsh.nlevels
            ndofs = ndofs + hspace.ndof_per_level(ilev);
            if (hmsh.nel_per_level(ilev) > 0)
                ind_e = (Ne(ilev)+1):Ne(ilev+1);
                sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
                b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, ...
                                reshape (aux(:, ind_e), [1, nqn_total, numel(ind_e)]));
                dofs = 1:ndofs;
                est(dofs) = est(dofs) + hspace.Csub{ilev}.' * b_lev;
            end
        end
        est = coef .* est;
        est = C0_est * sqrt (est);

    otherwise
        error('adaptivity_estimate_div_sigma_el:flag', ...
              'Unknown adaptivity_data.flag = %s', adaptivity_data.flag);
end

end


% =========================================================================
function divergence = accumulate_divergence(divergence, hgrad, ifield, ncomp)
% Accumulate the gradient of Voigt component ifield into the divergence
% vector.  hgrad has shape [ndim x nqn x nel].
%
% Voigt mapping:
%   1 -> sigma_xx: div(sigma)_1 += d(sigma_xx)/dx_1
%   2 -> sigma_yy: div(sigma)_2 += d(sigma_yy)/dx_2
%   3 -> sigma_zz: div(sigma)_3 += d(sigma_zz)/dx_3        (3D only)
%   4 -> sigma_xy: div(sigma)_1 += d(sigma_xy)/dx_2
%                  div(sigma)_2 += d(sigma_xy)/dx_1
%   5 -> sigma_yz: div(sigma)_2 += d(sigma_yz)/dx_3         (3D only)
%                  div(sigma)_3 += d(sigma_yz)/dx_2          (3D only)
%   6 -> sigma_xz: div(sigma)_1 += d(sigma_xz)/dx_3         (3D only)
%                  div(sigma)_3 += d(sigma_xz)/dx_1          (3D only)

    switch ifield
        case 1  % sigma_xx
            divergence(1,:,:) = divergence(1,:,:) + hgrad(1,:,:);
        case 2  % sigma_yy
            divergence(2,:,:) = divergence(2,:,:) + hgrad(2,:,:);
        case 3  % sigma_zz  (3D only)
            if ncomp == 3
                divergence(3,:,:) = divergence(3,:,:) + hgrad(3,:,:);
            end
        case 4  % sigma_xy (symmetric: sigma_xy = sigma_yx)
            divergence(1,:,:) = divergence(1,:,:) + hgrad(2,:,:);
            divergence(2,:,:) = divergence(2,:,:) + hgrad(1,:,:);
        case 5  % sigma_yz (3D only)
            if ncomp == 3
                divergence(2,:,:) = divergence(2,:,:) + hgrad(3,:,:);
                divergence(3,:,:) = divergence(3,:,:) + hgrad(2,:,:);
            end
        case 6  % sigma_xz (3D only)
            if ncomp == 3
                divergence(1,:,:) = divergence(1,:,:) + hgrad(3,:,:);
                divergence(3,:,:) = divergence(3,:,:) + hgrad(1,:,:);
            end
    end
end
