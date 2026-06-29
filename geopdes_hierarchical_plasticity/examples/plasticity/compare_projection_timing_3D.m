% COMPARE_PROJECTION_TIMING_3D  Wall-clock comparison of L2, QI, and QI_C0
% projection methods for the 3D rotating cylinder with centrifugal body force.
%
% Uses only the centrifugal body force (no pressure) to keep the problem
% purely elastic, allowing a clean timing comparison without Newton iteration
% variability.
%
% Output:
%   Console table of wall-clock times per projection method
%   results/projection_timing_3D.mat

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
% PHYSICAL PARAMETERS — centrifugal only (elastic, clean timing)
% =========================================================================
E  = 210000;   % [MPa]
nu = 0.3;
sigma_y = 240; % [MPa]
R_i = 100; R_o = 200; L = 50;  % [mm]

rho_omega2 = 5.0e-3;  % [MPa/mm^2]

fprintf('3D Projection Timing Comparison\n');
fprintf('================================\n');
fprintf('  Centrifugal body force only (rho*omega^2 = %.4e)\n', rho_omega2);
fprintf('  No internal pressure — problem stays elastic for clean timing\n\n');

%% ========================================================================
% PROBLEM DATA (3D)
% =========================================================================
nload = 4;  % fewer load steps for timing

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
method_data_template.nsub_coarse    = [8 8 2];
method_data_template.nsub_refine    = [2 2 2];
method_data_template.nquad          = [p+1 p+1 p+1];
method_data_template.space_type     = 'standard';
method_data_template.truncated      = 1;
method_data_template.nload          = nload;
method_data_template.newton_tol     = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;

adaptivity_data_template.flag = 'elements';
adaptivity_data_template.estimator = 'div_sigma';  % same estimator for all
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
% RUN EACH PROJECTION METHOD
% =========================================================================
proj_methods = {'L2', 'QI', 'QI_C0'};
T = struct();
T.methods = proj_methods;
T.nload = nload;
T.wall_total    = nan(1, 3);
T.ndof_final    = nan(1, 3);
T.nel_final     = nan(1, 3);
T.nlevels_final = nan(1, 3);

for imethod = 1:numel(proj_methods)
    method_name = proj_methods{imethod};
    fprintf('\n================================================================\n');
    fprintf('  Projection method: %s\n', method_name);
    fprintf('================================================================\n');

    method_data = method_data_template;
    method_data.type_projection = method_name;
    adaptivity_data = adaptivity_data_template;

    wall0 = tic;
    try
        [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
         cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
            adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);

        T.wall_total(imethod) = toc(wall0);
        T.ndof_final(imethod) = cell_hspace{end}.ndof;
        T.nel_final(imethod)  = cell_hmsh{end}.nel;
        T.nlevels_final(imethod) = cell_hmsh{end}.nlevels;

        fprintf('  Total wall time: %.1f s\n', T.wall_total(imethod));
        fprintf('  Final DOFs: %d, Elements: %d, Levels: %d\n', ...
                T.ndof_final(imethod), T.nel_final(imethod), T.nlevels_final(imethod));
    catch e
        T.wall_total(imethod) = toc(wall0);
        fprintf('  FAILED after %.1f s: %s\n', T.wall_total(imethod), e.message);
        T.ndof_final(imethod) = NaN;
        T.nel_final(imethod)  = NaN;
        T.nlevels_final(imethod) = NaN;
    end
end

%% ========================================================================
% RESULTS SUMMARY
% =========================================================================
fprintf('\n\n');
fprintf('=====================================================================\n');
fprintf('  3D PROJECTION TIMING COMPARISON  (%d load steps, div_sigma est.)\n', nload);
fprintf('=====================================================================\n');
fprintf('  Method    |  Wall [s]  |  DOFs   |  Elements  |  Levels\n');
fprintf('-------------------------------------------------------------\n');
for k = 1:3
    if isnan(T.ndof_final(k))
        fprintf('  %-8s  |  %8.1f  |  FAILED |            |\n', ...
                proj_methods{k}, T.wall_total(k));
    else
        fprintf('  %-8s  |  %8.1f  |  %5d  |  %5d     |  %d\n', ...
                proj_methods{k}, T.wall_total(k), ...
                T.ndof_final(k), T.nel_final(k), T.nlevels_final(k));
    end
end
fprintf('=====================================================================\n');

% Speedup ratios
if ~isnan(T.wall_total(2))
    fprintf('\n  QI / L2 wall-time ratio:     %.1fx slower\n', T.wall_total(2)/T.wall_total(1));
end
if ~isnan(T.wall_total(3))
    fprintf('  QI_C0 / L2 wall-time ratio:  %.1fx slower\n', T.wall_total(3)/T.wall_total(1));
end
if ~isnan(T.wall_total(2)) && ~isnan(T.wall_total(3))
    fprintf('  QI_C0 / QI wall-time ratio:  %.1fx\n', T.wall_total(3)/T.wall_total(2));
end

%% ========================================================================
% SAVE
% =========================================================================
save(fullfile(results_dir, 'projection_timing_3D.mat'), '-struct', 'T');

% CSV
fid = fopen(fullfile(results_dir, 'projection_timing_3D.csv'), 'w');
fprintf(fid, 'method,wall_s,ndof,nel,nlevels\n');
for k = 1:3
    fprintf(fid, '%s,%.3f,%d,%d,%d\n', ...
            proj_methods{k}, T.wall_total(k), T.ndof_final(k), T.nel_final(k), T.nlevels_final(k));
end
fclose(fid);

fprintf('\nResults saved to %s\n', results_dir);
