% DIAGNOSE_QI_SLOWNESS  Quick diagnostic: run one load step with QI and QI_C0
% and measure wall-clock + capture pinv fallback diagnostics from stderr.

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');

addpath(genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'initialize')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
addpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

cd(here);
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

%% Physical parameters
E  = 210000;   nu = 0.3;   sigma_y = 240;
rho_omega2 = 5.0e-3;

%% Problem data
problem_data.geo_name = fullfile(here, '..', 'data_files', 'geo_quarter_cylinder_3D.txt');
problem_data.nmnn_sides   = [1 2];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [];
problem_data.symm_sides   = [3 4 5 6];
problem_data.slider_sides = [];

mu_lame    = E / (2*(1+nu));
kappa_lame = E / (3*(1-2*nu));
problem_data.yield_stress = @(x, y, z) sigma_y * ones(size(x));
problem_data.kappa_lame   = @(x, y, z) kappa_lame * ones(size(x));
problem_data.mu_lame      = @(x, y, z) mu_lame * ones(size(x));

problem_data.f = @(x, y, z) rho_omega2 * cat(1, ...
    reshape(x, [1, size(x,1), size(x,2), size(x,3)]), ...
    reshape(y, [1, size(x,1), size(x,2), size(x,3)]), ...
    zeros([1, size(x,1), size(x,2), size(x,3)]));
problem_data.p = @(x, y, z) 0 * ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(3, size(x,1), size(x,2), size(x,3));
problem_data.h = @(x, y, z, ind, mult) zeros(3, size(x,1), size(x,2), size(x,3));

%% Test 1 load step, no adaptivity, on [4 4 1] mesh
p = 2;
method_data.degree         = [p p p];
method_data.regularity     = [p-1 p-1 p-1];
method_data.nsub_coarse    = [4 4 1];
method_data.nsub_refine    = [2 2 2];
method_data.nquad          = [p+1 p+1 p+1];
method_data.space_type     = 'standard';
method_data.truncated      = 1;
method_data.nload          = 1;
method_data.newton_tol     = 1e-8;
method_data.newton_tol_abs = 1e-10;
method_data.newton_iter_max = 100;

adaptivity_data.flag = 'elements';
adaptivity_data.estimator = 'div_sigma';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.5;
adaptivity_data.mark_param_coarsening = 0;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_ndof  = 200000;
adaptivity_data.max_nel   = 100000;
adaptivity_data.num_max_iter = 1;
adaptivity_data.tol = 1e-5 * 1e6;
adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any';
adaptivity_data.adm = p;
adaptivity_data.max_level = 1;

%% Run all three methods and report
for method_name = {'L2', 'QI', 'QI_C0'}
    mn = method_name{1};
    method_data.type_projection = mn;

    fprintf('\n=== %s projection ===\n', mn);
    wall0 = tic;
    [~, cell_hmsh, cell_hspace, cell_hspace_scalar, ~, ~, ~, ~] = ...
        adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);
    t = toc(wall0);

    fprintf('  Wall time: %.2f s\n', t);
    fprintf('  ndof=%d, nel=%d, scalar_ndof=%d\n', ...
        cell_hspace{end}.ndof, cell_hmsh{end}.nel, cell_hspace_scalar{end}.ndof);
end

%% Also test on [8 8 2] mesh
fprintf('\n\n======== LARGER MESH [8 8 2] ========\n');
method_data.nsub_coarse = [8 8 2];

for method_name = {'L2', 'QI', 'QI_C0'}
    mn = method_name{1};
    method_data.type_projection = mn;

    fprintf('\n=== %s projection ([8 8 2]) ===\n', mn);
    wall0 = tic;
    [~, cell_hmsh, cell_hspace, cell_hspace_scalar, ~, ~, ~, ~] = ...
        adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);
    t = toc(wall0);

    fprintf('  Wall time: %.2f s\n', t);
    fprintf('  ndof=%d, nel=%d, scalar_ndof=%d\n', ...
        cell_hspace{end}.ndof, cell_hmsh{end}.nel, cell_hspace_scalar{end}.ndof);
end

fprintf('\n\nDiagnostics above (from stderr) show pinv fallback counts.\n');
fprintf('Check MATLAB command window for qi_localLS_compute messages.\n');
