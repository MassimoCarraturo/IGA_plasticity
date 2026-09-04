% COMPARE_PROJECTION_ACCURACY_3D  Accuracy comparison of L2, QI, and QI_C0
% projection methods for the 3D rotating cylinder with centrifugal body force.
%
% Uses only the centrifugal body force (no pressure) to keep the problem
% purely elastic.  Compares each adaptive method against a fine uniform-mesh
% reference solution.
%
% Metrics:
%   - Radial displacement error at outer boundary (each load step)
%   - Radial stress profile error along radial line (theta=45 deg, mid-height)
%   - Hoop stress profile error along radial line
%
% Output:
%   Console table + results/projection_accuracy_3D.mat + .csv

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
results_dir = fullfile(here, 'results');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
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

%% ========================================================================
% PHYSICAL PARAMETERS — centrifugal only (elastic)
% =========================================================================
E  = 210000;   % [MPa]
nu = 0.3;
sigma_y = 240; % [MPa]
R_i = 100; R_o = 200; L = 50;  % [mm]

rho_omega2 = 5.0e-3;  % [MPa/mm^2]

fprintf('3D Projection Accuracy Comparison\n');
fprintf('==================================\n');
fprintf('  Centrifugal body force only (rho*omega^2 = %.4e)\n', rho_omega2);
fprintf('  No internal pressure — problem stays elastic\n\n');

%% ========================================================================
% PROBLEM DATA (3D)
% =========================================================================
nload = 4;

problem_data.geo_name = fullfile(here, '..', 'data_files', 'geo_quarter_cylinder_3D.txt');
problem_data.nmnn_sides   = [1 2];        % both inner and outer traction-free
problem_data.drchlt_sides = [];
problem_data.press_sides  = [];           % no pressure
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

problem_data.p = @(x, y, z) 0 * ones(size(x));  % no pressure
problem_data.g = @(x, y, z, ind) zeros(3, size(x,1), size(x,2), size(x,3));
problem_data.h = @(x, y, z, ind, mult) zeros(3, size(x,1), size(x,2), size(x,3));

%% ========================================================================
% DISCRETIZATION — common template (3D)
% =========================================================================
p = 2;
method_data_template.degree         = [p p p];
method_data_template.regularity     = [p-1 p-1 p-1];
method_data_template.nsub_coarse    = [4 4 1];
method_data_template.nsub_refine    = [2 2 2];
method_data_template.nquad          = [p+1 p+1 p+1];
method_data_template.space_type     = 'standard';
method_data_template.truncated      = 1;
method_data_template.nload          = nload;
method_data_template.newton_tol     = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;

adaptivity_data_template.flag = 'elements';
adaptivity_data_template.estimator = 'div_sigma';
adaptivity_data_template.C0_est = 1.0;
adaptivity_data_template.mark_param = 0.5;
adaptivity_data_template.mark_param_coarsening = 0.01;
adaptivity_data_template.mark_strategy = 'MS';
adaptivity_data_template.max_ndof  = 200000;
adaptivity_data_template.max_nel   = 100000;
adaptivity_data_template.num_max_iter = 3;
adaptivity_data_template.tol = 1e-5 * 1e6;
adaptivity_data_template.adm_strategy = 'admissible';
adaptivity_data_template.coarsening_flag = 'any';
adaptivity_data_template.adm = p;
adaptivity_data_template.max_level = 3;

%% ========================================================================
% STEP 1: LOAD REFERENCE SOLUTION (pre-computed by compute_reference_3D.m)
% =========================================================================
ref_file = fullfile(results_dir, 'reference_3D.mat');
if ~exist(ref_file, 'file')
    error('Reference file not found: %s\nRun compute_reference_3D.m first.', ref_file);
end
ref = load(ref_file);
u_r_ref_vec = ref.u_r_ref;
sigma_r_ref = ref.sigma_r_ref;
sigma_t_ref = ref.sigma_t_ref;
rad_pos     = ref.rad_pos;
n_eval      = ref.n_eval;

fprintf('  Loaded reference: ndof=%d, nel=%d (computed in %.1fs)\n', ...
        ref.ref_ndof, ref.ref_nel, ref.ref_time);
fprintf('  Reference stress range: sigma_r = [%.2f, %.2f], sigma_t = [%.2f, %.2f] MPa\n', ...
        min(sigma_r_ref(:,end)), max(sigma_r_ref(:,end)), ...
        min(sigma_t_ref(:,end)), max(sigma_t_ref(:,end)));

%% ========================================================================
% STEP 2: ADAPTIVE SOLVES — L2, QI, QI_C0
% =========================================================================
proj_methods = {'L2', 'QI', 'QI_C0'};
n_methods = numel(proj_methods);

R = struct();
R.methods = proj_methods;
R.nload = nload;
R.n_eval = n_eval;
R.radius = linspace(R_i, R_o, n_eval);

fields = {'wall','ndof','nel','nlevels', ...
           'err_u_max','err_u_avg', ...
           'err_sigma_r_max','err_sigma_r_avg', ...
           'err_sigma_t_max','err_sigma_t_avg'};
for k = 1:numel(fields)
    R.(fields{k}) = nan(1, n_methods);
end
R.sigma_r = cell(1, n_methods);
R.sigma_t = cell(1, n_methods);
R.u_r     = cell(1, n_methods);

for im = 1:n_methods
    method_name = proj_methods{im};
    fprintf('\n================================================================\n');
    fprintf('  ADAPTIVE: projection = %s\n', method_name);
    fprintf('================================================================\n');

    method_data = method_data_template;
    method_data.type_projection = method_name;
    adaptivity_data = adaptivity_data_template;

    wall0 = tic;
    try
        [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
         cell_u, ~, cell_sigma, sol_data] = ...
            adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);
    catch e
        fprintf('  FAILED: %s\n', e.message);
        R.wall(im) = toc(wall0);
        continue;
    end
    R.wall(im) = toc(wall0);
    R.ndof(im)    = cell_hspace{end}.ndof;
    R.nel(im)     = cell_hmsh{end}.nel;
    R.nlevels(im) = cell_hmsh{end}.nlevels;

    % Evaluate displacement at outer boundary
    u_r_im = zeros(1, nload);
    for i = 1:nload
        [eu, ~] = sp_eval(cell_u{i}, cell_hspace{i}, geometry, {1, 0.5, 0.5});
        u_r_im(i) = norm(eu(1:2));
    end
    R.u_r{im} = u_r_im;

    % Evaluate stress along radial line
    sigma_r_im = zeros(n_eval, nload);
    sigma_t_im = zeros(n_eval, nload);
    for ls = 1:nload
        for jpt = 1:n_eval
            position = {rad_pos(jpt), 0.5, 0.5};
            sigma_xx = sp_eval(cell_sigma{ls}(:,1), cell_hspace_scalar{ls}, geometry, position);
            sigma_yy = sp_eval(cell_sigma{ls}(:,2), cell_hspace_scalar{ls}, geometry, position);
            sigma_xy = sp_eval(cell_sigma{ls}(:,4), cell_hspace_scalar{ls}, geometry, position);
            c = sqrt(2)/2;
            sigma_r_im(jpt, ls) = c^2*(sigma_xx + sigma_yy) + 2*c^2*sigma_xy;
            sigma_t_im(jpt, ls) = c^2*(sigma_xx + sigma_yy) - 2*c^2*sigma_xy;
        end
    end
    R.sigma_r{im} = sigma_r_im;
    R.sigma_t{im} = sigma_t_im;

    % Compute errors vs reference
    diff_u = abs(u_r_im - u_r_ref_vec);
    R.err_u_max(im) = max(diff_u);
    R.err_u_avg(im) = mean(diff_u);

    diff_sr = abs(sigma_r_im - sigma_r_ref);
    diff_st = abs(sigma_t_im - sigma_t_ref);
    R.err_sigma_r_max(im) = max(diff_sr(:));
    R.err_sigma_r_avg(im) = mean(diff_sr(:));
    R.err_sigma_t_max(im) = max(diff_st(:));
    R.err_sigma_t_avg(im) = mean(diff_st(:));

    fprintf('  ndof=%d  nel=%d  nlev=%d  wall=%.1fs\n', ...
            R.ndof(im), R.nel(im), R.nlevels(im), R.wall(im));
    fprintf('  err_u   max=%.4e  avg=%.4e\n', R.err_u_max(im), R.err_u_avg(im));
    fprintf('  sigma_r max=%.4e  avg=%.4e\n', R.err_sigma_r_max(im), R.err_sigma_r_avg(im));
    fprintf('  sigma_t max=%.4e  avg=%.4e\n', R.err_sigma_t_max(im), R.err_sigma_t_avg(im));
end

%% ========================================================================
% STEP 3: RESULTS SUMMARY
% =========================================================================
fprintf('\n\n');
fprintf('========================================================================\n');
fprintf('  3D PROJECTION ACCURACY COMPARISON  (%d load steps, div_sigma est.)\n', nload);
fprintf('========================================================================\n');
fprintf('  Metric              |  %-14s |  %-14s |  %-14s\n', proj_methods{:});
fprintf('------------------------------------------------------------------------\n');
fprintf('  DOFs                |  %12d  |  %12d  |  %12d\n', R.ndof(1), R.ndof(2), R.ndof(3));
fprintf('  Elements            |  %12d  |  %12d  |  %12d\n', R.nel(1), R.nel(2), R.nel(3));
fprintf('  Levels              |  %12d  |  %12d  |  %12d\n', R.nlevels(1), R.nlevels(2), R.nlevels(3));
fprintf('  Wall time [s]       |  %12.1f  |  %12.1f  |  %12.1f\n', R.wall(1), R.wall(2), R.wall(3));
fprintf('  err_u max           |  %12.4e  |  %12.4e  |  %12.4e\n', R.err_u_max(1), R.err_u_max(2), R.err_u_max(3));
fprintf('  err_u avg           |  %12.4e  |  %12.4e  |  %12.4e\n', R.err_u_avg(1), R.err_u_avg(2), R.err_u_avg(3));
fprintf('  err_sigma_r max     |  %12.4e  |  %12.4e  |  %12.4e\n', R.err_sigma_r_max(1), R.err_sigma_r_max(2), R.err_sigma_r_max(3));
fprintf('  err_sigma_r avg     |  %12.4e  |  %12.4e  |  %12.4e\n', R.err_sigma_r_avg(1), R.err_sigma_r_avg(2), R.err_sigma_r_avg(3));
fprintf('  err_sigma_t max     |  %12.4e  |  %12.4e  |  %12.4e\n', R.err_sigma_t_max(1), R.err_sigma_t_max(2), R.err_sigma_t_max(3));
fprintf('  err_sigma_t avg     |  %12.4e  |  %12.4e  |  %12.4e\n', R.err_sigma_t_avg(1), R.err_sigma_t_avg(2), R.err_sigma_t_avg(3));
fprintf('========================================================================\n');

% Per-metric winner
metrics = {'err_u_max','err_u_avg','err_sigma_r_max','err_sigma_r_avg','err_sigma_t_max','err_sigma_t_avg'};
wins = zeros(1, n_methods);
for k = 1:numel(metrics)
    vals = R.(metrics{k});
    [~, idx] = min(vals);
    wins(idx) = wins(idx) + 1;
end
fprintf('\n');
for im = 1:n_methods
    fprintf('  %s wins %d / %d metrics\n', proj_methods{im}, wins(im), numel(metrics));
end

% Efficiency: error per DOF
fprintf('\n  Efficiency (err_sigma_r_avg * ndof):\n');
for im = 1:n_methods
    if ~isnan(R.ndof(im))
        fprintf('    %s: %.2e\n', proj_methods{im}, R.err_sigma_r_avg(im) * R.ndof(im));
    end
end

%% ========================================================================
% SAVE
% =========================================================================
R.sigma_r_ref = sigma_r_ref;
R.sigma_t_ref = sigma_t_ref;
R.u_r_ref = u_r_ref_vec;
R.ref_ndof = ref.ref_ndof;
R.ref_nel  = ref.ref_nel;
R.rad_pos  = rad_pos;

save(fullfile(results_dir, 'projection_accuracy_3D.mat'), '-struct', 'R');

fid = fopen(fullfile(results_dir, 'projection_accuracy_3D.csv'), 'w');
fprintf(fid, 'method,wall_s,ndof,nel,nlevels,err_u_max,err_u_avg,err_sr_max,err_sr_avg,err_st_max,err_st_avg\n');
for im = 1:n_methods
    fprintf(fid, '%s,%.3f,%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
            proj_methods{im}, R.wall(im), R.ndof(im), R.nel(im), R.nlevels(im), ...
            R.err_u_max(im), R.err_u_avg(im), ...
            R.err_sigma_r_max(im), R.err_sigma_r_avg(im), ...
            R.err_sigma_t_max(im), R.err_sigma_t_avg(im));
end
fclose(fid);

fprintf('\nResults saved to %s\n', results_dir);
