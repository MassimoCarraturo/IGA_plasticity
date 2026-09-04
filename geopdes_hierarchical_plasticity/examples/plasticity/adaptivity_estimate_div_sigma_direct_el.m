function est = adaptivity_estimate_div_sigma_direct_el (u, eps_pl_store, ...
    sigma_store, geometry, hmsh, hspace, hspace_scalar, ...
    problem_data, adaptivity_data)
% ADAPTIVITY_ESTIMATE_DIV_SIGMA_DIRECT_EL  Element-local residual error
% indicator based on the strong-form equilibrium residual, computed
% DIRECTLY from the displacement u and the projected plastic strain eps_pl.
%
% Volume term:
%   eta_K^2 = h_K^2 * || f + div(sigma) ||^2_{L^2(K)}
%
% where div(sigma) is computed analytically via the constitutive law:
%   sigma = lambda * tr(eps - eps_pl) I + 2*mu * (eps - eps_pl)
%   div(sigma) = (lambda+mu)*grad(div u) + mu*laplacian(u)
%                - lambda*grad(tr eps_pl) - 2*mu*div(eps_pl)
%
% This avoids projecting sigma onto a C^1 space: second derivatives of u
% (C^1 displacement) and first derivatives of eps_pl (any continuity) are
% evaluated element-wise for L^2 integration.
%
% Boundary term:
%   + sum_{F in dK cap Gamma_N} h_F * || g - sigma_h . n ||^2_{L^2(F)}
%
% The boundary traction uses sigma_store (pointwise evaluation — no
% differentiation, so any continuity works).
%
% Voigt convention: 1=xx, 2=yy, 3=zz, 4=xy, 5=yz, 6=xz
%
% INPUT:
%   u              - displacement DOFs [hspace.ndof x 1]
%   eps_pl_store   - projected plastic strain DOFs [hspace_scalar.ndof x 6]
%   sigma_store    - projected stress DOFs [hspace_scalar.ndof x 6]
%                    (used for Neumann boundary term only)
%   geometry       - NURBS geometry structure
%   hmsh           - hierarchical mesh (primal / displacement)
%   hspace         - hierarchical displacement space (vector, C^1)
%   hspace_scalar  - hierarchical scalar space (for eps_pl)
%   problem_data   - physical data (f, g, mu_lame, kappa_lame, ...)
%   adaptivity_data- estimator parameters (C0_est, flag, ...)
%
% OUTPUT:
%   est            - element-wise error indicators

if (isfield(adaptivity_data, 'C0_est'))
    C0_est = adaptivity_data.C0_est;
else
    C0_est = 1;
end

ncomp      = hspace.ncomp;          % 2 or 3
ndim       = hmsh.ndim;
nqn_total  = hmsh.mesh_of_level(1).nqn;

% =========================================================================
%  VOLUME TERM: div(sigma) from second derivatives of u and gradients of eps_pl
% =========================================================================
divergence = zeros(ncomp, nqn_total, hmsh.nel);
F = zeros(hmsh.rdim, nqn_total, hmsh.nel);

last_dof_u = cumsum(hspace.ndof_per_level);
last_dof_s = cumsum(hspace_scalar.ndof_per_level);

el_shift = 0;
for ilev = 1:hmsh.nlevels
    nel_lev = hmsh.nel_per_level(ilev);
    if (nel_lev == 0); continue; end

    msh_lev = hmsh.msh_lev{ilev};
    el_range = el_shift + (1:nel_lev);

    % --- Evaluate displacement Hessian: H[comp, dir1, dir2, qn, el] -----
    sp_u_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, ...
        'value', false, 'gradient', false, 'hessian', true);
    u_lev = hspace.Csub{ilev} * u(1:last_dof_u(ilev));
    H_u = sp_eval_msh (u_lev, sp_u_lev, msh_lev, 'hessian');
    % H_u: [ncomp, ndim, ndim, nqn, nel_lev]

    % --- Evaluate eps_pl gradients: grad_epl[voigt, dir, qn, el] --------
    sp_s_lev = sp_evaluate_element_list (hspace_scalar.space_of_level(ilev), msh_lev, ...
        'value', false, 'gradient', true);
    grad_epl = zeros(6, ndim, nqn_total, nel_lev);
    for c = 1:6
        epl_lev = hspace_scalar.Csub{ilev} * eps_pl_store(1:last_dof_s(ilev), c);
        tmp = sp_eval_msh (epl_lev, sp_s_lev, msh_lev, 'gradient');
        % tmp: [ndim, nqn, nel_lev]
        grad_epl(c, :, :, :) = tmp;
    end

    % --- Material properties at Gauss points -----------------------------
    x = cell(hmsh.rdim, 1);
    for idim = 1:hmsh.rdim
        x{idim} = reshape(msh_lev.geo_map(idim,:,:), msh_lev.nqn, nel_lev);
    end
    mu_val    = reshape(problem_data.mu_lame(x{:}),    [1, nqn_total, nel_lev]);
    kappa_val = reshape(problem_data.kappa_lame(x{:}),  [1, nqn_total, nel_lev]);
    lambda_val = kappa_val - 2*mu_val/3;

    % --- Assemble div(sigma) at quad points of this level ----------------
    % div(sigma)_i = (lambda+mu)*grad_div_u_i + mu*laplacian_u_i
    %              - lambda*grad_tr_epl_i - 2*mu*div_epl_i
    for i = 1:ncomp
        % (lambda+mu) * d(div u)/dx_i = (lambda+mu) * sum_j d^2 u_j / (dx_j dx_i)
        grad_div_u_i = zeros(1, nqn_total, nel_lev);
        for j = 1:ndim
            grad_div_u_i = grad_div_u_i + reshape(H_u(j, j, i, :, :), [1, nqn_total, nel_lev]);
        end

        % mu * laplacian(u_i) = mu * sum_j d^2 u_i / dx_j^2
        laplacian_u_i = zeros(1, nqn_total, nel_lev);
        for j = 1:ndim
            laplacian_u_i = laplacian_u_i + reshape(H_u(i, j, j, :, :), [1, nqn_total, nel_lev]);
        end

        % lambda * d(tr eps_pl)/dx_i   (Voigt 1+2+3 = xx+yy+zz)
        grad_tr_epl_i = reshape(grad_epl(1,i,:,:) + grad_epl(2,i,:,:), [1, nqn_total, nel_lev]);
        if ndim == 3
            grad_tr_epl_i = grad_tr_epl_i + reshape(grad_epl(3,i,:,:), [1, nqn_total, nel_lev]);
        end

        % 2*mu * sum_j d(eps_pl_{ij})/dx_j   (divergence of eps_pl row i)
        div_epl_i = compute_epl_div_row(grad_epl, i, ndim, nqn_total, nel_lev);

        divergence(i, :, el_range) = ...
              (lambda_val + mu_val) .* grad_div_u_i ...
            + mu_val .* laplacian_u_i ...
            - lambda_val .* grad_tr_epl_i ...
            - 2*mu_val .* div_epl_i;
    end

    F(:, :, el_range) = msh_lev.geo_map;
    el_shift = el_shift + nel_lev;
end

% =========================================================================
%  BODY FORCE
% =========================================================================
x = cell(hmsh.rdim, 1);
for idim = 1:hmsh.rdim
    x{idim} = reshape(F(idim, :, :), [], hmsh.nel);
end
valf = problem_data.f(x{:});

% Residual: (f + div sigma)^2, summed over components
aux = sum((valf + divergence).^2, 1);
aux = reshape(aux, nqn_total, hmsh.nel);

% =========================================================================
%  ELEMENT-WISE INDICATOR
% =========================================================================
switch adaptivity_data.flag
    case 'elements'
        w = []; h = [];
        for ilev = 1:hmsh.nlevels
            if (hmsh.msh_lev{ilev}.nel ~= 0)
                w = cat(2, w, hmsh.msh_lev{ilev}.quad_weights .* hmsh.msh_lev{ilev}.jacdet);
                h = cat(1, h, hmsh.msh_lev{ilev}.element_size(:));
            end
        end
        h = h * sqrt(hmsh.ndim);

        % Integrate ||f + div sigma||^2 over each element
        elem_int = sum(aux .* w, 1);
        est_sq = (h(:)).^2 .* elem_int(:);

        % --- Neumann boundary term ---
        bnd_contrib = compute_neumann_boundary_term(sigma_store, ...
            hmsh, hmsh, hspace_scalar, problem_data, ncomp, size(sigma_store,2));
        est_sq = est_sq + bnd_contrib;

        est = C0_est * sqrt(est_sq);

    otherwise
        error('adaptivity_estimate_div_sigma_direct_el:flag', ...
              'Only flag=''elements'' is supported. Got: %s', adaptivity_data.flag);
end

end


% =========================================================================
function div_epl_i = compute_epl_div_row(grad_epl, i, ndim, nqn, nel)
% Divergence of the eps_pl tensor, row i:
%   sum_j  d(eps_pl_{ij}) / dx_j
%
% Voigt mapping (3D): 1=xx 2=yy 3=zz 4=xy 5=yz 6=xz
%   eps_pl tensor row 1: [xx xy xz] = Voigt [1 4 6]
%   eps_pl tensor row 2: [xy yy yz] = Voigt [4 2 5]
%   eps_pl tensor row 3: [xz yz zz] = Voigt [6 5 3]
%
% Voigt mapping (2D): 1=xx 2=yy 3=xy
%   eps_pl tensor row 1: [xx xy] = Voigt [1 3]
%   eps_pl tensor row 2: [xy yy] = Voigt [3 2]

if ndim == 3
    voigt_row = {[1 4 6], [4 2 5], [6 5 3]};
elseif ndim == 2
    voigt_row = {[1 3], [3 2]};
end

row = voigt_row{i};
div_epl_i = zeros(1, nqn, nel);
for j = 1:ndim
    div_epl_i = div_epl_i + reshape(grad_epl(row(j), j, :, :), [1, nqn, nel]);
end

end


% =========================================================================
function bnd_contrib = compute_neumann_boundary_term (sigma_store, ...
    hmsh, hmsh_scalar, hspace_scalar, problem_data, ncomp, ncomp_sigma)
% Neumann boundary contribution (reused from adaptivity_estimate_div_sigma_el).
%
%   bnd_contrib(K) = sum_{F in dK cap Gamma_N}  h_F * int_F |g - sigma_h . n|^2 dF
%
% Uses sigma_store for pointwise evaluation — no differentiation needed,
% so any continuity (including C^0) works.

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

    % Evaluate boundary geometry (normals, weights, jacdet)
    msh_side = hmsh_eval_boundary_side(hmsh_scalar, iside);
    nqn_bnd = msh_side.nqn;
    nel_bnd = msh_side.nel;

    % Evaluate all Voigt stress components at boundary quad pts
    sigma_at_bnd = zeros(ncomp_sigma, nqn_bnd, nel_bnd);
    bnd_shift = cumsum([0; hmsh_bnd.nel_per_level(:)]);

    for ilev = 1:hmsh_scalar.nlevels
        if (ilev > numel(hmsh_bnd.nel_per_level) || hmsh_bnd.nel_per_level(ilev) == 0)
            continue
        end

        hmsh_scalar.mesh_of_level(ilev).boundary = ...
            hmsh_scalar.boundary(1).mesh_of_level(ilev);
        for is = 2:numel(hmsh_scalar.boundary)
            hmsh_scalar.mesh_of_level(ilev).boundary(is) = ...
                hmsh_scalar.boundary(iside).mesh_of_level(ilev);
        end

        msh_side_from_int = msh_boundary_side_from_interior( ...
            hmsh_scalar.mesh_of_level(ilev), iside);
        sp_bnd = hspace_scalar.space_of_level(ilev).constructor(msh_side_from_int);
        msh_int_struct = msh_evaluate_element_list(msh_side_from_int, ...
            hmsh_bnd.active{ilev});
        sp_bnd = sp_evaluate_element_list(sp_bnd, msh_int_struct, 'value', true);

        bnd_range = (bnd_shift(ilev)+1):bnd_shift(ilev+1);

        for ifield = 1:ncomp_sigma
            u_lev = hspace_scalar.Csub{ilev} * ...
                sigma_store(1:last_dof_scalar(ilev), ifield);
            sigma_val = sp_eval_msh(u_lev, sp_bnd, msh_int_struct, 'value');
            sigma_at_bnd(ifield, :, bnd_range) = reshape(sigma_val, 1, nqn_bnd, []);
        end
    end

    % Compute traction t_i = sum_j sigma_ij * n_j
    normals  = msh_side.normal;
    traction = voigt_traction(sigma_at_bnd, normals, ncomp);

    % Prescribed Neumann traction g
    x_bnd = cell(hmsh.rdim, 1);
    for idim = 1:hmsh.rdim
        x_bnd{idim} = reshape(msh_side.geo_map(idim,:,:), nqn_bnd, nel_bnd);
    end
    if isfield(problem_data, 'g')
        g_val = problem_data.g(x_bnd{:}, iside);
    else
        g_val = zeros(ncomp, nqn_bnd, nel_bnd);
    end

    % Integrate |g - t|^2 over each boundary face
    resid = sum((g_val - traction).^2, 1);
    resid = reshape(resid, nqn_bnd, nel_bnd);

    bnd_weights  = msh_side.quad_weights .* msh_side.jacdet;
    bnd_integral = sum(resid .* bnd_weights, 1)';

    h_F = msh_side.element_size(:) * sqrt(ndim);

    % Map boundary elements -> volume elements and accumulate
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
[~, nqn, nel] = size(normals);
traction = zeros(ncomp, nqn, nel);

if ncomp == 3
    traction(1,:,:) = sigma_at_bnd(1,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(4,:,:) .* normals(2,:,:) ...
                    + sigma_at_bnd(6,:,:) .* normals(3,:,:);
    traction(2,:,:) = sigma_at_bnd(4,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(2,:,:) .* normals(2,:,:) ...
                    + sigma_at_bnd(5,:,:) .* normals(3,:,:);
    traction(3,:,:) = sigma_at_bnd(6,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(5,:,:) .* normals(2,:,:) ...
                    + sigma_at_bnd(3,:,:) .* normals(3,:,:);
elseif ncomp == 2
    traction(1,:,:) = sigma_at_bnd(1,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(3,:,:) .* normals(2,:,:);
    traction(2,:,:) = sigma_at_bnd(3,:,:) .* normals(1,:,:) ...
                    + sigma_at_bnd(2,:,:) .* normals(2,:,:);
end
end


% =========================================================================
function bnd_to_vol = map_bnd_to_vol_elements (hmsh, hmsh_scalar, iside)
hmsh_bnd = hmsh_scalar.boundary(iside);
ndim     = hmsh_scalar.ndim;
ind2     = ceil(iside/2);
ind      = setdiff(1:ndim, ind2);

bnd_to_vol = zeros(hmsh_bnd.nel, 1);
bnd_shift  = cumsum([0; hmsh_bnd.nel_per_level(:)]);
vol_shift  = cumsum([0; hmsh.nel_per_level(:)]);

for ilev = 1:hmsh_bnd.nlevels
    if hmsh_bnd.nel_per_level(ilev) == 0, continue; end

    nel_dir     = hmsh_scalar.mesh_of_level(ilev).nel_dir;
    bnd_nel_dir = nel_dir(ind);
    bnd_active  = hmsh_bnd.active{ilev};
    vol_active  = hmsh.active{ilev};

    if mod(iside, 2) == 1
        fixed_sub = 1;
    else
        fixed_sub = nel_dir(ind2);
    end

    for i = 1:numel(bnd_active)
        bnd_subs = cell(ndim-1, 1);
        [bnd_subs{:}] = ind2sub([bnd_nel_dir, 1], bnd_active(i));

        vol_subs = cell(ndim, 1);
        for d = 1:ndim-1
            vol_subs{ind(d)} = bnd_subs{d};
        end
        vol_subs{ind2} = fixed_sub;

        vol_el = sub2ind([nel_dir, 1], vol_subs{:});
        pos = find(vol_active == vol_el, 1);
        if ~isempty(pos)
            bnd_to_vol(bnd_shift(ilev) + i) = vol_shift(ilev) + pos;
        end
    end
end
end
