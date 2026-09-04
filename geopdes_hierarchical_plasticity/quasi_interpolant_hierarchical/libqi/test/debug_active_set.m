% DEBUG_ACTIVE_SET  Compare sp_get_basis_functions with the C-equivalent
% (TP product over [a, b+deg]) for every active basis function.

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..', '..');

addpath(genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));

deg  = [3 3];
nelx = 8; nely = 8;
geometry = geo_load(nrb4surf([0 0], [1 0], [0 1], [1 1]));
[knots, zeta] = kntrefine(geometry.nurbs.knots, [nelx nely] - 1, deg, deg - 1);
rule     = msh_gauss_nodes(deg + 1);
[qn, qw] = msh_set_quad_nodes(zeta, rule);
msh   = msh_cartesian(zeta, qn, qw, geometry);
space = sp_bspline(knots, deg, msh);
hmsh   = hierarchical_mesh(msh, [2 2]);
hspace = hierarchical_space(hmsh, space, 'standard', false, deg - 1);

sp_lv  = hspace.space_of_level(1);
ms_lv  = hmsh.mesh_of_level(1);
nd     = sp_lv.ndof_dir;
nel    = ms_lv.nel_dir;
fprintf('ndof_dir = [%d %d], nel_dir = [%d %d]\n', nd(1), nd(2), nel(1), nel(2));

n_diff = 0;
for k = 1:hspace.ndof
    kl = hspace.active{1}(k);
    [kx, ky] = ind2sub(nd, kl);
    a_x = max(1, kx - deg(1)); b_x = min(nel(1), kx);
    a_y = max(1, ky - deg(2)); b_y = min(nel(2), ky);
    support = sp_get_cells(sp_lv, ms_lv, kl);
    ind_matlab = sort(sp_get_basis_functions(sp_lv, ms_lv, support));

    rng_x = a_x:min(nd(1), b_x + deg(1));
    rng_y = a_y:min(nd(2), b_y + deg(2));
    [I, J] = ndgrid(rng_x, rng_y);
    ind_mine = sort(I(:) + (J(:) - 1) * nd(1));

    if ~isequal(ind_matlab(:), ind_mine(:))
        n_diff = n_diff + 1;
        if n_diff <= 5
            fprintf('k=%d kl=%d: matlab len=%d, mine len=%d\n', k, kl, ...
                    numel(ind_matlab), numel(ind_mine));
            fprintf('  matlab: %s\n', mat2str(ind_matlab(:)'));
            fprintf('  mine:   %s\n', mat2str(ind_mine(:)'));
        end
    end
end
fprintf('Total active-set mismatches: %d / %d\n', n_diff, hspace.ndof);

% Now also check support cells
n_diff_supp = 0;
for k = 1:hspace.ndof
    kl = hspace.active{1}(k);
    [kx, ky] = ind2sub(nd, kl);
    a_x = max(1, kx - deg(1)); b_x = min(nel(1), kx);
    a_y = max(1, ky - deg(2)); b_y = min(nel(2), ky);
    support_matlab = sort(sp_get_cells(sp_lv, ms_lv, kl));
    [I, J] = ndgrid(a_x:b_x, a_y:b_y);
    support_mine = sort(I(:) + (J(:) - 1) * nel(1));
    if ~isequal(support_matlab(:), support_mine(:))
        n_diff_supp = n_diff_supp + 1;
        if n_diff_supp <= 5
            fprintf('SUPPORT k=%d: matlab=%s, mine=%s\n', k, ...
                    mat2str(support_matlab(:)'), mat2str(support_mine(:)'));
        end
    end
end
fprintf('Total support mismatches: %d / %d\n', n_diff_supp, hspace.ndof);
