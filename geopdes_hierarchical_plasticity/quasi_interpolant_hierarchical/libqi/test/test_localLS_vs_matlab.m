% TEST_LOCALLS_VS_MATLAB  Correctness test: getcoeff_localLS_Bspl_c vs the
% reference getcoeff_localLS_Bspl on a single-level 2D scattered-data fit.

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..', '..');

addpath(genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
addpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

% Build a tiny single-level 2D space + mesh
deg  = [3 3];
nelx = 8; nely = 8;
geometry = geo_load (nrb4surf([0 0], [1 0], [0 1], [1 1]));
[knots, zeta] = kntrefine(geometry.nurbs.knots, [nelx nely] - 1, deg, deg - 1);
rule     = msh_gauss_nodes(deg + 1);
[qn, qw] = msh_set_quad_nodes(zeta, rule);
msh   = msh_cartesian(zeta, qn, qw, geometry);
space = sp_bspline(knots, deg, msh);
hmsh   = hierarchical_mesh(msh, [2 2]);
hspace = hierarchical_space(hmsh, space, 'standard', false, deg - 1);

% Regularity = 0 (Bezier-like, used by QI_ref) — separate hspace2
[knots0, zeta0] = kntrefine(geometry.nurbs.knots, [nelx nely] - 1, deg, [0 0]);
[qn0, qw0] = msh_set_quad_nodes(zeta0, rule);
msh0   = msh_cartesian(zeta0, qn0, qw0, geometry);
space0 = sp_bspline(knots0, deg, msh0);
hmsh0   = hierarchical_mesh(msh0, [2 2]);
hspace0 = hierarchical_space(hmsh0, space0, 'standard', false, [0 0]);
fprintf('hspace (default reg) ndof = %d, hspace0 (reg=0) ndof = %d\n', hspace.ndof, hspace0.ndof);

% Quasi-uniform "quadrature" data: many points per element so that the
% local LS systems are well-conditioned (mirrors the actual use case in
% history_variable_projection_hier where data comes from Gauss nodes).
rng(0);
[gx, gy] = ndgrid(linspace(0.01, 0.99, 50), linspace(0.01, 0.99, 50));
x = gx(:) + 0.005 * (rand(numel(gx), 1) - 0.5);
y = gy(:) + 0.005 * (rand(numel(gy), 1) - 0.5);
M = numel(x);
ncomp = 3;
fv = [(tanh(9*y - 9*x) + 1)/9, exp(-((10*x-6).^2 + (10*y+7).^2))/1.5, sin(3*pi*x).*cos(2*pi*y)];
data = [x, y];
fprintf('M = %d data points\n', M);

warning ('off', 'MATLAB:nearlySingularMatrix');
warning ('off', 'MATLAB:singularMatrix');

for lambda = [0, 1e-9, 1e-3]
    fprintf('--- lambda = %g ---\n', lambda);

    tic;
    QI_ref = getcoeff_localLS_Bspl(hspace, hmsh, data, fv, lambda);
    t_ref = toc;

    tic;
    QI_c = getcoeff_localLS_Bspl_c(hspace, hmsh, data, fv, lambda);
    t_c = toc;

    err = max(abs(QI_ref(:) - QI_c(:)));
    rel = err / max(abs(QI_ref(:)));
    fprintf('  ndof=%d  ncomp=%d   max|ref - c| = %.3e  (rel %.3e)\n', ...
            hspace.ndof, ncomp, err, rel);
    fprintf('  reference time = %.3fs   libqi time = %.3fs   speedup x%.2f\n', ...
            t_ref, t_c, t_ref / max(t_c, eps));

    if rel < 1e-8
        fprintf('  PASS\n');
    else
        fprintf('  FAIL\n');
    end
end

% Same comparison on the regularity-0 space (mirrors the QI_ref projection)
fprintf('\n========= regularity-0 space (QI_ref-like) =========\n');
for lambda = [0, 1e-9, 1e-3]
    fprintf('--- lambda = %g ---\n', lambda);
    QI_ref0 = getcoeff_localLS_Bspl(hspace0, hmsh0, data, fv, lambda);
    QI_c0   = getcoeff_localLS_Bspl_c(hspace0, hmsh0, data, fv, lambda);
    err = max(abs(QI_ref0(:) - QI_c0(:)));
    rel = err / max(abs(QI_ref0(:)));
    fprintf('  ndof=%d  max|ref - c| = %.3e  (rel %.3e)\n', hspace0.ndof, err, rel);
    if rel < 1e-8, fprintf('  PASS\n'); else, fprintf('  FAIL\n'); end
end

disp('Test done.');
