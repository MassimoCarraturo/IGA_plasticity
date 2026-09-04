% TEST_NEUMANN_BOUNDARY_TERM  Detailed diagnostics for the Neumann boundary
% term in the div_sigma estimator.
%
% Runs the solver once, then calls the estimator twice:
%   (a) with nmnn_sides to get volume + boundary
%   (b) with empty nmnn_sides to get volume-only
% and reports the per-element breakdown.

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

%% Mesh sizes to test
mesh_configs = {[4 4 1], [8 8 2]};

for im = 1:numel(mesh_configs)
    nsub = mesh_configs{im};
    fprintf('\n================================================================\n');
    fprintf('  Mesh [%d %d %d]\n', nsub);
    fprintf('================================================================\n');

    p = 2;
    method_data.degree         = [p p p];
    method_data.regularity     = [p-1 p-1 p-1];
    method_data.nsub_coarse    = nsub;
    method_data.nsub_refine    = [2 2 2];
    method_data.nquad          = [p+1 p+1 p+1];
    method_data.space_type     = 'standard';
    method_data.truncated      = 1;
    method_data.nload          = 1;
    method_data.newton_tol     = 1e-8;
    method_data.newton_tol_abs = 1e-10;
    method_data.newton_iter_max = 100;
    method_data.type_projection = 'L2';

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

    % Run the solver
    [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
     ~, ~, cell_sigma, ~] = ...
        adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);

    hmsh = cell_hmsh{end};
    hspace = cell_hspace{end};
    hspace_scalar = cell_hspace_scalar{end};
    sigma_store = cell_sigma{end};  % [ndof_scalar x 6] Voigt stress DOFs
    hmsh_scalar = hmsh;  % same when num_bisections = 0

    % (a) Full estimator (volume + boundary)
    est_full = adaptivity_estimate_div_sigma_el( ...
        sigma_store, geometry, hmsh, hspace, hmsh_scalar, hspace_scalar, ...
        problem_data, adaptivity_data);

    % (b) Volume-only estimator (remove nmnn_sides)
    pd_novol = problem_data;
    pd_novol.nmnn_sides = [];
    est_vol = adaptivity_estimate_div_sigma_el( ...
        sigma_store, geometry, hmsh, hspace, hmsh_scalar, hspace_scalar, ...
        pd_novol, adaptivity_data);

    % Compute boundary-only contribution
    est_bnd_sq = est_full.^2 - est_vol.^2;

    % Report
    eta_vol  = sqrt(sum(est_vol.^2));
    eta_full = sqrt(sum(est_full.^2));
    eta_bnd  = sqrt(sum(max(est_bnd_sq, 0)));
    pct_bnd  = 100 * eta_bnd^2 / eta_full^2;

    fprintf('\n  Total eta (volume only): %.6e\n', eta_vol);
    fprintf('  Total eta (vol+bnd):     %.6e\n', eta_full);
    fprintf('  Boundary contribution:   %.6e  (%.1f%% of total squared)\n', eta_bnd, pct_bnd);

    n_with_bnd = sum(est_bnd_sq > 0);
    fprintf('  Elements with boundary term: %d / %d\n', n_with_bnd, hmsh.nel);

    % Show elements with largest boundary contribution
    [sorted_bnd, idx] = sort(est_bnd_sq, 'descend');
    n_show = min(8, n_with_bnd);
    if n_show > 0
        fprintf('\n  Top %d elements by boundary contribution:\n', n_show);
        fprintf('  %6s  %12s  %12s  %12s  %8s\n', ...
            'elem', 'eta_vol', 'eta_bnd', 'eta_full', 'bnd/tot%');
        for i = 1:n_show
            ie = idx(i);
            bnd_pct = 100 * max(est_bnd_sq(ie),0) / est_full(ie)^2;
            fprintf('  %6d  %12.4e  %12.4e  %12.4e  %7.1f%%\n', ...
                ie, est_vol(ie), sqrt(max(est_bnd_sq(ie),0)), est_full(ie), bnd_pct);
        end
    end
end

fprintf('\n\nAll tests completed successfully.\n');
