function est = adaptivity_estimate_div_sigma_el (sigma_store, geometry, hmsh, hspace, hmsh_scalar, hspace_scalar, problem_data, adaptivity_data)
% ADAPTIVITY_ESTIMATE_DIV_SIGMA_EL  Element-local residual error indicator
% based on the strong-form equilibrium residual:
%
%   eta_K^2 = h_K^2 * || f + div(sigma_h) ||^2_{L^2(K)}
%           + sum_{F in dK cap Gamma_N} h_F * || g - sigma_h . n ||^2_{L^2(F)}
%
% The first (volume) term measures how well the discrete stress locally
% satisfies equilibrium.  The second (boundary) term penalises the mismatch
% between the computed traction and the prescribed Neumann data g on every
% boundary face F that belongs to a Neumann side.
%
% Voigt convention: 1=sig_xx, 2=sig_yy, 3=sig_zz, 4=sig_xy, 5=sig_yz, 6=sig_xz
%
% Two evaluation paths for the volume term:
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
        est_sq = (h(:)).^2 .* elem_int(:);

        % --- Neumann boundary term: h_F * || g - sigma_h . n ||^2 ---
        bnd_contrib = compute_neumann_boundary_term(sigma_store, ...
            hmsh, hmsh_scalar, hspace_scalar, problem_data, ncomp, ncomp_sigma);
        est_sq = est_sq + bnd_contrib;

        est = C0_est * sqrt (est_sq);

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


% =========================================================================
function bnd_contrib = compute_neumann_boundary_term (sigma_store, ...
    hmsh, hmsh_scalar, hspace_scalar, problem_data, ncomp, ncomp_sigma)
% Compute the Neumann boundary contribution to the element-wise estimator:
%
%   bnd_contrib(K) = sum_{F in dK cap Gamma_N}  h_F * int_F |g - sigma_h . n|^2 dF
%
% The sum runs over Neumann boundary faces F of volume element K.
% Returns a vector of size [hmsh.nel x 1].

bnd_contrib = zeros(hmsh.nel, 1);

if ~isfield(problem_data, 'nmnn_sides') || isempty(problem_data.nmnn_sides)
    return
end

ndim = hmsh.ndim;
last_dof_scalar = cumsum(hspace_scalar.ndof_per_level);

for iside = problem_data.nmnn_sides
    if iside > numel(hmsh_scalar.boundary), continue; end
    hmsh_bnd = hmsh_scalar.boundary(iside);
    if hmsh_bnd.nel == 0, continue; end

    % ----- Evaluate boundary geometry (normals, weights, jacdet) ---------
    msh_side = hmsh_eval_boundary_side(hmsh_scalar, iside);
    nqn_bnd = msh_side.nqn;
    nel_bnd = msh_side.nel;

    % ----- Evaluate all Voigt stress components at boundary quad pts -----
    sigma_at_bnd = zeros(ncomp_sigma, nqn_bnd, nel_bnd);
    bnd_shift = cumsum([0; hmsh_bnd.nel_per_level(:)]);

    for ilev = 1:hmsh_scalar.nlevels
        if (ilev > numel(hmsh_bnd.nel_per_level) || hmsh_bnd.nel_per_level(ilev) == 0)
            continue
        end

        % Populate mesh_of_level.boundary so msh_eval_boundary_side works
        % (boundary was cleared in the constructor for memory; we restore
        %  the entry for side iside here — local copy, does not persist)
        hmsh_scalar.mesh_of_level(ilev).boundary = ...
            hmsh_scalar.boundary(1).mesh_of_level(ilev);
        for is = 2:numel(hmsh_scalar.boundary)
            hmsh_scalar.mesh_of_level(ilev).boundary(is) = ...
                hmsh_scalar.boundary(iside).mesh_of_level(ilev);
        end

        % Interior-view boundary mesh and space
        msh_side_from_int = msh_boundary_side_from_interior( ...
            hmsh_scalar.mesh_of_level(ilev), iside);
        sp_bnd = hspace_scalar.space_of_level(ilev).constructor(msh_side_from_int);
        msh_int_struct = msh_evaluate_element_list(msh_side_from_int, ...
            hmsh_bnd.active{ilev});
        sp_bnd = sp_evaluate_element_list(sp_bnd, msh_int_struct, 'value', true);

        % Boundary-element range in the global (across levels) list
        bnd_range = (bnd_shift(ilev)+1):bnd_shift(ilev+1);

        for ifield = 1:ncomp_sigma
            u_lev = hspace_scalar.Csub{ilev} * ...
                sigma_store(1:last_dof_scalar(ilev), ifield);
            sigma_val = sp_eval_msh(u_lev, sp_bnd, msh_int_struct, 'value');
            % sigma_val: [1 x nqn_bnd x nel_bnd_lev] for scalar space
            sigma_at_bnd(ifield, :, bnd_range) = reshape(sigma_val, 1, nqn_bnd, []);
        end
    end

    % ----- Compute traction t_i = sum_j sigma_ij * n_j ------------------
    normals  = msh_side.normal;                         % [rdim x nqn x nel]
    traction = voigt_traction(sigma_at_bnd, normals, ncomp);  % [ncomp x nqn x nel]

    % ----- Prescribed Neumann traction g ---------------------------------
    x_bnd = cell(hmsh.rdim, 1);
    for idim = 1:hmsh.rdim
        x_bnd{idim} = reshape(msh_side.geo_map(idim,:,:), nqn_bnd, nel_bnd);
    end
    if isfield(problem_data, 'g')
        g_val = problem_data.g(x_bnd{:}, iside);       % [ncomp x nqn x nel]
    else
        g_val = zeros(ncomp, nqn_bnd, nel_bnd);
    end

    % ----- Integrate |g - t|^2 over each boundary face -------------------
    resid = sum((g_val - traction).^2, 1);              % [1 x nqn x nel]
    resid = reshape(resid, nqn_bnd, nel_bnd);

    bnd_weights  = msh_side.quad_weights .* msh_side.jacdet;  % [nqn x nel]
    bnd_integral = sum(resid .* bnd_weights, 1)';       % [nel_bnd x 1]

    % Face element size (scaled consistently with volume h)
    h_F = msh_side.element_size(:) * sqrt(ndim);        % [nel_bnd x 1]

    % ----- Map boundary elements -> volume elements and accumulate -------
    bnd_to_vol = map_bnd_to_vol_elements(hmsh, hmsh_scalar, iside);

    for ib = 1:nel_bnd
        iv = bnd_to_vol(ib);
        if iv > 0
            bnd_contrib(iv) = bnd_contrib(iv) + h_F(ib) * bnd_integral(ib);
        end
    end
end

end


% =========================================================================
function traction = voigt_traction (sigma_at_bnd, normals, ncomp)
% Compute traction  t_i = sum_j sigma_ij * n_j  from Voigt stress.
%
%   sigma_at_bnd : [ncomp_sigma x nqn x nel]   Voigt stress at boundary pts
%   normals      : [rdim x nqn x nel]           unit outward normals
%   ncomp        : 2 (2D) or 3 (3D)
%
% Voigt order (3D): 1=xx 2=yy 3=zz 4=xy 5=yz 6=xz
% Voigt order (2D): 1=xx 2=yy 3=xy

[~, nqn, nel] = size(normals);
traction = zeros(ncomp, nqn, nel);

if ncomp == 3
    % t_1 = sig_xx*n1 + sig_xy*n2 + sig_xz*n3
    traction(1,:,:) = sigma_at_bnd(1,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(4,:,:) .* normals(2,:,:) ...
                    + sigma_at_bnd(6,:,:) .* normals(3,:,:);
    % t_2 = sig_xy*n1 + sig_yy*n2 + sig_yz*n3
    traction(2,:,:) = sigma_at_bnd(4,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(2,:,:) .* normals(2,:,:) ...
                    + sigma_at_bnd(5,:,:) .* normals(3,:,:);
    % t_3 = sig_xz*n1 + sig_yz*n2 + sig_zz*n3
    traction(3,:,:) = sigma_at_bnd(6,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(5,:,:) .* normals(2,:,:) ...
                    + sigma_at_bnd(3,:,:) .* normals(3,:,:);

elseif ncomp == 2
    % t_1 = sig_xx*n1 + sig_xy*n2
    traction(1,:,:) = sigma_at_bnd(1,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(3,:,:) .* normals(2,:,:);
    % t_2 = sig_xy*n1 + sig_yy*n2
    traction(2,:,:) = sigma_at_bnd(3,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(2,:,:) .* normals(2,:,:);
end

end


% =========================================================================
function bnd_to_vol = map_bnd_to_vol_elements (hmsh, hmsh_scalar, iside)
% Map global boundary element indices to global volume element indices.
%
%   bnd_to_vol(ib) = global volume-element index for boundary element ib,
%                    or 0 if the parent volume element is not active
%                    (should not happen in a consistent mesh).
%
% Works level-by-level: at each level, converts the boundary element index
% (in the (d-1)-dimensional boundary mesh) to the adjacent volume element
% index (in the d-dimensional mesh) using the structured tensor-product layout.

hmsh_bnd = hmsh_scalar.boundary(iside);
ndim     = hmsh_scalar.ndim;
ind2     = ceil(iside/2);                      % fixed parametric direction
ind      = setdiff(1:ndim, ind2);              % free directions

bnd_to_vol = zeros(hmsh_bnd.nel, 1);
bnd_shift  = cumsum([0; hmsh_bnd.nel_per_level(:)]);
vol_shift  = cumsum([0; hmsh.nel_per_level(:)]);

for ilev = 1:hmsh_bnd.nlevels
    if hmsh_bnd.nel_per_level(ilev) == 0, continue; end

    nel_dir     = hmsh_scalar.mesh_of_level(ilev).nel_dir;
    bnd_nel_dir = nel_dir(ind);     % boundary mesh dimensions
    bnd_active  = hmsh_bnd.active{ilev};
    vol_active  = hmsh.active{ilev};

    % Value of the fixed subscript: 1 (low side) or max (high side)
    if mod(iside, 2) == 1
        fixed_sub = 1;
    else
        fixed_sub = nel_dir(ind2);
    end

    for i = 1:numel(bnd_active)
        % Boundary element -> subscripts in boundary mesh
        bnd_subs = cell(ndim-1, 1);
        [bnd_subs{:}] = ind2sub([bnd_nel_dir, 1], bnd_active(i));

        % Assemble volume subscripts
        vol_subs = cell(ndim, 1);
        for d = 1:ndim-1
            vol_subs{ind(d)} = bnd_subs{d};
        end
        vol_subs{ind2} = fixed_sub;

        vol_el = sub2ind([nel_dir, 1], vol_subs{:});

        % Find in active volume elements at this level
        pos = find(vol_active == vol_el, 1);
        if ~isempty(pos)
            bnd_to_vol(bnd_shift(ilev) + i) = vol_shift(ilev) + pos;
        end
    end
end

end
