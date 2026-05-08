% SMOKE_TEST_SPHERE_QI  Reduced-mesh QI run to verify libqi backend works
% end-to-end inside history_variable_projection_hier.

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

nload = 2;                                % small for smoke test
P     = 332*.99;
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.p = @(x, y, z) P*ones (size (x));

p = 2;
method_data.degree         = [p p p];
method_data.regularity     = [p-1 p-1 p-1];
method_data.nsub_coarse    = [10, 4, 4];
method_data.nsub_refine    = [2 2 2];
method_data.nquad          = [p+1 p+1 p+1];
method_data.space_type     = 'standard';
method_data.truncated      = 1;
method_data.nload          = nload;
method_data.newton_tol     = 1e-6;
method_data.newton_tol_abs = 1e-8;
method_data.newton_iter_max = 50;
method_data.type_projection = 'QI';

adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.9;
adaptivity_data.mark_param_coarsening = 0.05;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 1;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 1;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5 * 3.14 * 30000;
adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any';
adaptivity_data.adm = p;

if exist('qi_local_ls_mex', 'file') == 3
    fprintf('libqi qi_local_ls_mex is on path: backend will be used.\n');
else
    fprintf('libqi qi_local_ls_mex NOT found: falling back to MATLAB reference.\n');
end

tic;
[geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
 cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
    adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);
fprintf('SMOKE TEST QI OK in %.2f s\n', toc);
fprintf('|cell_eps_pl{end}|_max = %.4e\n', max(abs(cell_eps_pl{end}(:))));
fprintf('|cell_sigma{end}|_max  = %.4e\n', max(abs(cell_sigma{end}(:))));
