% SMOKE_TEST_ESTIMATOR  Quick run of QI / max_level=2 with both
% estimators to verify the new stress-gradient indicator works and to
% compare which elements each one targets.

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
addpath (genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath (genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'initialize')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
addpath (fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

cd (here);
warning ('off', 'MATLAB:nearlySingularMatrix');
warning ('off', 'MATLAB:singularMatrix');

problem_data.geo_name = 'geo_eighth_sphere.txt';
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [1];
problem_data.symm_sides   = [3 4 5];
problem_data.slider_sides = [6];
problem_data.penalty_slider = @(x, y, z) 1e8 * ones (size (x));
E  = 210000; nu = 0.3;
problem_data.yield_stress = @(x, y, z) 240 * ones (size (x));
problem_data.kappa_lame   = @(x, y, z) E/(3*(1-2*nu)) * ones (size (x));
problem_data.mu_lame      = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));
nload = 3; P = 332*.99;
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.p = @(x, y, z) P*ones (size (x));

p = 2;
method_data.degree         = [p p p];
method_data.regularity     = [p-1 p-1 p-1];
method_data.nsub_coarse    = [3, 1, 1];
method_data.nsub_refine    = [2 2 2];
method_data.nquad          = [p+1 p+1 p+1];
method_data.space_type     = 'standard';
method_data.truncated      = 1;
method_data.nload          = nload;
method_data.newton_tol     = 1e-7;
method_data.newton_tol_abs = 1e-9;
method_data.newton_iter_max = 50;
method_data.type_projection = 'QI';

adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.5;
adaptivity_data.mark_param_coarsening = 0.0;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 2;
adaptivity_data.max_ndof = 200000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 100000;
adaptivity_data.tol = 1e-5 * 3.14 * 30000;
adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any';
adaptivity_data.adm = p;

for est_kind = {'div_sigma', 'stress_gradient'}
    fprintf('\n==================== %s ====================\n', est_kind{1});
    adaptivity_data.estimator = est_kind{1};

    cpu0 = cputime(); wall0 = tic;
    [~, cell_hmsh, ~, ~, ~, ~, ~, ~] = ...
        adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);
    fprintf('  wall = %.2f s   CPU = %.2f s\n', toc(wall0), cputime() - cpu0);
    fprintf('  final mesh: %d elements at %d levels\n', cell_hmsh{end}.nel, cell_hmsh{end}.nlevels);
end
