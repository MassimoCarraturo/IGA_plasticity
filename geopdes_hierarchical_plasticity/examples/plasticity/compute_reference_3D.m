% COMPUTE_REFERENCE_3D  Compute and save the fine-mesh reference solution.
%
% Saves: results/reference_3D.mat

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

%% Physical parameters
E  = 210000;   nu = 0.3;   sigma_y = 240;
R_i = 100; R_o = 200; L = 50;
rho_omega2 = 5.0e-3;
nload = 4;

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

%% Discretization
p = 2;
method_data_ref.degree         = [p p p];
method_data_ref.regularity     = [p-1 p-1 p-1];
method_data_ref.nsub_coarse    = [32 32 8];
method_data_ref.nsub_refine    = [2 2 2];
method_data_ref.nquad          = [p+1 p+1 p+1];
method_data_ref.space_type     = 'standard';
method_data_ref.truncated      = 1;
method_data_ref.nload          = nload;
method_data_ref.newton_tol     = 1e-8;
method_data_ref.newton_tol_abs = 1e-10;
method_data_ref.newton_iter_max = 100;
method_data_ref.type_projection = 'L2';

adaptivity_data_ref.flag = 'elements';
adaptivity_data_ref.estimator = 'div_sigma';
adaptivity_data_ref.C0_est = 1.0;
adaptivity_data_ref.mark_param = 0.5;
adaptivity_data_ref.mark_param_coarsening = 0;
adaptivity_data_ref.mark_strategy = 'MS';
adaptivity_data_ref.max_ndof  = 200000;
adaptivity_data_ref.max_nel   = 100000;
adaptivity_data_ref.num_max_iter = 1;
adaptivity_data_ref.tol = 1e-5 * 1e6;
adaptivity_data_ref.adm_strategy = 'admissible';
adaptivity_data_ref.coarsening_flag = 'any';
adaptivity_data_ref.adm = p;
adaptivity_data_ref.max_level = 1;

%% Solve reference
fprintf('Computing reference solution: [32 32 8], %d load steps ...\n', nload);
wall0 = tic;
[geometry_ref, cell_hmsh_ref, cell_hspace_ref, cell_hspace_scalar_ref, ...
 cell_u_ref, ~, cell_sigma_ref, ~] = ...
    adaptivity_J2_plasticity(problem_data, method_data_ref, adaptivity_data_ref);
ref_time = toc(wall0);
fprintf('  Reference: ndof=%d, nel=%d, time=%.1fs\n', ...
        cell_hspace_ref{end}.ndof, cell_hmsh_ref{end}.nel, ref_time);

%% Evaluate reference quantities
n_eval = 100;
rad_pos = linspace(0, 1, n_eval);

% Displacement at outer boundary
u_r_ref_vec = zeros(1, nload);
for i = 1:nload
    [eu, ~] = sp_eval(cell_u_ref{i}, cell_hspace_ref{i}, geometry_ref, {1, 0.5, 0.5});
    u_r_ref_vec(i) = norm(eu(1:2));
end

% Stress along radial line (theta=45 deg, mid-height)
sigma_r_ref = zeros(n_eval, nload);
sigma_t_ref = zeros(n_eval, nload);
for ls = 1:nload
    for jpt = 1:n_eval
        position = {rad_pos(jpt), 0.5, 0.5};
        sigma_xx = sp_eval(cell_sigma_ref{ls}(:,1), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        sigma_yy = sp_eval(cell_sigma_ref{ls}(:,2), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        sigma_xy = sp_eval(cell_sigma_ref{ls}(:,4), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        c = sqrt(2)/2;
        sigma_r_ref(jpt, ls) = c^2*(sigma_xx + sigma_yy) + 2*c^2*sigma_xy;
        sigma_t_ref(jpt, ls) = c^2*(sigma_xx + sigma_yy) - 2*c^2*sigma_xy;
    end
end

fprintf('  sigma_r range: [%.2f, %.2f], sigma_t range: [%.2f, %.2f] MPa\n', ...
        min(sigma_r_ref(:,end)), max(sigma_r_ref(:,end)), ...
        min(sigma_t_ref(:,end)), max(sigma_t_ref(:,end)));

%% Save
ref = struct();
ref.u_r_ref     = u_r_ref_vec;
ref.sigma_r_ref = sigma_r_ref;
ref.sigma_t_ref = sigma_t_ref;
ref.rad_pos     = rad_pos;
ref.n_eval      = n_eval;
ref.nload       = nload;
ref.ref_ndof    = cell_hspace_ref{end}.ndof;
ref.ref_nel     = cell_hmsh_ref{end}.nel;
ref.ref_time    = ref_time;
save(fullfile(results_dir, 'reference_3D.mat'), '-struct', 'ref');
fprintf('Reference saved to %s\n', fullfile(results_dir, 'reference_3D.mat'));
