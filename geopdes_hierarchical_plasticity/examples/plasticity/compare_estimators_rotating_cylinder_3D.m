% COMPARE_ESTIMATORS_ROTATING_CYLINDER_3D  3D version of the estimator
% comparison for a thick-walled rotating cylinder with internal pressure.
%
% Problem: 3D quarter-cylinder (R_i=100, R_o=200, L=50 mm) in generalized
% plane strain (z-symmetry on both axial faces), loaded by:
%   (1) Internal pressure P at the inner radius
%   (2) Centrifugal body force f = rho*omega^2 * [x; y; 0] (z-axis rotation)
%
% The same plane-strain analytical solution applies (validated against 2D).
% This script demonstrates that the div_sigma estimator advantage
% persists in 3D, where the equilibrium residual ||f + div(sigma)||
% now involves a 3D divergence operator and 6 stress components.
%
% Boundary conditions (quarter-cylinder with z-symmetry):
%   Side 1 (r=R_i):    internal pressure
%   Side 2 (r=R_o):    traction-free
%   Sides 3, 4:        circumferential symmetry (x-z and y-z planes)
%   Sides 5, 6 (z=0,L): axial symmetry (uz=0, plane strain)
%
% Reference solution: fine uniform mesh (no adaptivity).

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
% PHYSICAL PARAMETERS
% =========================================================================
E  = 210000;   % Young's modulus [MPa]
nu = 0.3;      % Poisson's ratio
sigma_y = 240; % Yield stress [MPa]
R_i = 100;     % Inner radius [mm]
R_o = 200;     % Outer radius [mm]
L   = 50;      % Cylinder height [mm] (arbitrary, plane-strain in z)

% --- Combined loading: internal pressure + centrifugal body force ---
P = 85;           % [MPa] internal pressure
rho_omega2 = 5.0e-3;  % [MPa/mm^2] centrifugal body force coefficient

% Verify elastic stress levels (same as 2D — plane strain)
[sig_r_el, sig_t_el, ~, r_ref] = rotating_cylinder_solution(E, nu, rho_omega2, R_i, R_o, 200);
sigma_eq_centrifugal = sqrt(3)/2 * abs(sig_t_el - sig_r_el);
sigma_t_P = P * (R_o^2 + R_i^2) / (R_o^2 - R_i^2);
sigma_eq_pressure = sqrt(3)/2 * (sigma_t_P + P);

fprintf('3D Rotating Cylinder — Estimator Comparison\n');
fprintf('============================================\n');
fprintf('  Geometry: R_i=%g, R_o=%g, L=%g mm (quarter, plane-strain in z)\n', R_i, R_o, L);
fprintf('  P = %.1f MPa,  rho*omega^2 = %.4e MPa/mm^2\n', P, rho_omega2);
fprintf('  Centrifugal sigma_vM at R_i: %.1f MPa (%.0f%% of sigma_y)\n', ...
        sigma_eq_centrifugal(1), 100*sigma_eq_centrifugal(1)/sigma_y);
fprintf('  Pressure sigma_vM at R_i:    %.1f MPa (%.0f%% of sigma_y)\n', ...
        sigma_eq_pressure, 100*sigma_eq_pressure/sigma_y);
fprintf('  Combined (approx):           %.1f MPa (%.0f%% of sigma_y)\n\n', ...
        sigma_eq_centrifugal(1) + sigma_eq_pressure, ...
        100*(sigma_eq_centrifugal(1) + sigma_eq_pressure)/sigma_y);

%% ========================================================================
% PROBLEM DATA (3D)
% =========================================================================
nload = 8;  % load steps

problem_data.geo_name = fullfile(here, '..', 'data_files', 'geo_quarter_cylinder_3D.txt');

% Boundary sides:
%   1=inner (xi1=0), 2=outer (xi1=1)
%   3=theta=0 (xi2=0), 4=theta=90 (xi2=1)
%   5=z=0 (xi3=0), 6=z=L (xi3=1)
problem_data.nmnn_sides   = [2];           % outer: traction-free
problem_data.drchlt_sides = [];            % no Dirichlet
problem_data.press_sides  = [1];           % inner: pressure
problem_data.symm_sides   = [3 4 5 6];    % quarter symmetry + z-plane strain
problem_data.slider_sides = [];

% Material (3D function handles)
mu_lame    = E / (2*(1+nu));
kappa_lame = E / (3*(1-2*nu));
problem_data.yield_stress = @(x, y, z) sigma_y * ones(size(x));
problem_data.kappa_lame   = @(x, y, z) kappa_lame * ones(size(x));
problem_data.mu_lame      = @(x, y, z) mu_lame * ones(size(x));

% Body force: centrifugal f = rho*omega^2 * [x; y; 0] (rotation about z)
problem_data.f = @(x, y, z) rho_omega2 * cat(1, ...
    reshape(x, [1, size(x,1), size(x,2), size(x,3)]), ...
    reshape(y, [1, size(x,1), size(x,2), size(x,3)]), ...
    zeros([1, size(x,1), size(x,2), size(x,3)]));

% Pressure on inner surface
problem_data.p = @(x, y, z) P * ones(size(x));

% Traction-free outer boundary
problem_data.g = @(x, y, z, ind) zeros(3, size(x,1), size(x,2), size(x,3));
problem_data.h = @(x, y, z, ind, mult) zeros(3, size(x,1), size(x,2), size(x,3));

%% ========================================================================
% DISCRETIZATION — common template (3D)
% =========================================================================
p = 2;
method_data_template.degree         = [p p p];
method_data_template.regularity     = [p-1 p-1 p-1];
method_data_template.nsub_coarse    = [2 2 1];    % 2x2 in-plane, 1 in z (plane strain)
method_data_template.nsub_refine    = [2 2 2];
method_data_template.nquad          = [p+1 p+1 p+1];
method_data_template.space_type     = 'standard';
method_data_template.truncated      = 1;
method_data_template.nload          = nload;
method_data_template.newton_tol     = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;
method_data_template.type_projection = 'L2';  % L2 is much faster than QI in 3D

adaptivity_data_template.flag = 'elements';
adaptivity_data_template.C0_est = 1.0;
adaptivity_data_template.mark_param = 0.5;
adaptivity_data_template.mark_param_coarsening = 0.01;  % mild coarsening
adaptivity_data_template.mark_strategy = 'MS';
adaptivity_data_template.max_ndof  = 15000;   % allow more 3D refinement
adaptivity_data_template.max_nel   = 5000;
adaptivity_data_template.num_max_iter = 3;    % 3 adaptive iterations per load step
adaptivity_data_template.tol = 1e-5 * 1e6;   % effectively never reached
adaptivity_data_template.adm_strategy = 'admissible';
adaptivity_data_template.coarsening_flag = 'any';
adaptivity_data_template.adm = p;
adaptivity_data_template.max_level = 4;       % allow deeper refinement in 3D

%% ========================================================================
% STEP 1: REFERENCE SOLUTION (fine uniform mesh, no adaptivity)
% =========================================================================
fprintf('\n================================================================\n');
fprintf('  REFERENCE SOLUTION (3D fine uniform mesh)\n');
fprintf('================================================================\n');

method_data_ref = method_data_template;
method_data_ref.nsub_coarse = [8 8 1];        % fine in-plane, single layer in z
method_data_ref.type_projection = 'L2';

adaptivity_data_ref = adaptivity_data_template;
adaptivity_data_ref.num_max_iter = 1;         % no adaptivity
adaptivity_data_ref.max_level = 1;
adaptivity_data_ref.mark_param_coarsening = 0;

wall0 = tic;
[geometry_ref, cell_hmsh_ref, cell_hspace_ref, cell_hspace_scalar_ref, ...
 cell_u_ref, ~, cell_sigma_ref, ~] = ...
    adaptivity_J2_plasticity(problem_data, method_data_ref, adaptivity_data_ref);
ref_time = toc(wall0);
fprintf('  Reference: ndof=%d, nel=%d, time=%.1fs\n', ...
        cell_hspace_ref{end}.ndof, cell_hmsh_ref{end}.nel, ref_time);

% Evaluate reference stress along radial line (at z=0.5, theta=45 deg)
n_eval = 100;
rad_pos = linspace(0, 1, n_eval);
pt_eval_3D = {rad_pos, 0.5, 0.5};  % xi1=radial, xi2=0.5 (45 deg), xi3=0.5 (mid-height)

voigt_3D = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];
load_steps_eval = [round(nload/2), nload];

sigma_r_ref = zeros(n_eval, numel(load_steps_eval));
sigma_t_ref = zeros(n_eval, numel(load_steps_eval));
u_r_ref_vec = zeros(1, nload+1);

% Reference displacement at outer boundary (r=R_o, theta=45, z=0.5)
for i = 1:nload
    [eu, ~] = sp_eval(cell_u_ref{i}, cell_hspace_ref{i}, geometry_ref, {1, 0.5, 0.5});
    u_r_ref_vec(i+1) = norm(eu(1:2));  % radial displacement (x,y components)
end

% Reference stress along radial line
for jload = 1:numel(load_steps_eval)
    ls = load_steps_eval(jload);
    for jpt = 1:n_eval
        position = {rad_pos(jpt), 0.5, 0.5};
        % Evaluate all stress components
        sigma_xx = sp_eval(cell_sigma_ref{ls}(:,1), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        sigma_yy = sp_eval(cell_sigma_ref{ls}(:,2), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        sigma_xy = sp_eval(cell_sigma_ref{ls}(:,4), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        % At theta=45 deg: radial = [cos45, sin45], tangential = [-sin45, cos45]
        c = sqrt(2)/2;
        sigma_r_ref(jpt, jload) = c^2*(sigma_xx + sigma_yy) + 2*c^2*sigma_xy;
        sigma_t_ref(jpt, jload) = c^2*(sigma_xx + sigma_yy) - 2*c^2*sigma_xy;
    end
end

%% ========================================================================
% STEP 2: ADAPTIVE SOLVES — div_sigma vs stress_gradient
% =========================================================================
estimators = {'div_sigma', 'stress_gradient'};
proj_method = 'L2';

R = struct();
R.estimators = estimators;
R.proj_method = proj_method;
R.rho_omega2 = rho_omega2;
R.load_steps_eval = load_steps_eval;
R.radius = linspace(R_i, R_o, n_eval);
R.dimension = 3;

fields = {'wall','ndof','nel','nlevels', ...
           'err_u_max','err_u_avg', ...
           'err_sigma_r_max','err_sigma_r_avg', ...
           'err_sigma_t_max','err_sigma_t_avg'};
for k = 1:numel(fields)
    R.(fields{k}) = nan(1, numel(estimators));
end
R.sigma_r = cell(1, numel(estimators));
R.sigma_t = cell(1, numel(estimators));
R.u_r     = cell(1, numel(estimators));

for ie = 1:numel(estimators)
    fprintf('\n================================================================\n');
    fprintf('  3D ADAPTIVE: estimator = %s,  method = %s\n', estimators{ie}, proj_method);
    fprintf('================================================================\n');

    method_data = method_data_template;
    method_data.type_projection = proj_method;
    adaptivity_data = adaptivity_data_template;
    adaptivity_data.estimator = estimators{ie};

    wall0 = tic;
    try
        [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, cell_u, ~, cell_sigma, sol_data] = ...
            adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);
    catch e
        fprintf('  FAILED: %s\n', e.message);
        R.wall(ie) = toc(wall0);
        continue;
    end
    R.wall(ie) = toc(wall0);
    R.ndof(ie) = cell_hspace{end}.ndof;
    R.nel(ie)  = cell_hmsh{end}.nel;
    R.nlevels(ie) = cell_hmsh{end}.nlevels;

    % Evaluate displacement at outer boundary
    u_r_ie = zeros(1, nload+1);
    for i = 1:nload
        [eu, ~] = sp_eval(cell_u{i}, cell_hspace{i}, geometry, {1, 0.5, 0.5});
        u_r_ie(i+1) = norm(eu(1:2));
    end
    R.u_r{ie} = u_r_ie;

    % Evaluate stress along radial line
    sigma_r_ie = zeros(n_eval, numel(load_steps_eval));
    sigma_t_ie = zeros(n_eval, numel(load_steps_eval));
    for jload = 1:numel(load_steps_eval)
        ls = load_steps_eval(jload);
        for jpt = 1:n_eval
            position = {rad_pos(jpt), 0.5, 0.5};
            sigma_xx = sp_eval(cell_sigma{ls}(:,1), cell_hspace_scalar{ls}, geometry, position);
            sigma_yy = sp_eval(cell_sigma{ls}(:,2), cell_hspace_scalar{ls}, geometry, position);
            sigma_xy = sp_eval(cell_sigma{ls}(:,4), cell_hspace_scalar{ls}, geometry, position);
            c = sqrt(2)/2;
            sigma_r_ie(jpt, jload) = c^2*(sigma_xx + sigma_yy) + 2*c^2*sigma_xy;
            sigma_t_ie(jpt, jload) = c^2*(sigma_xx + sigma_yy) - 2*c^2*sigma_xy;
        end
    end
    R.sigma_r{ie} = sigma_r_ie;
    R.sigma_t{ie} = sigma_t_ie;

    % Compute errors vs reference
    diff_u = abs(u_r_ie(2:end) - u_r_ref_vec(2:end));
    R.err_u_max(ie) = max(diff_u);
    R.err_u_avg(ie) = mean(diff_u);

    diff_sr = abs(sigma_r_ie - sigma_r_ref);
    diff_st = abs(sigma_t_ie - sigma_t_ref);
    R.err_sigma_r_max(ie) = max(diff_sr(:));
    R.err_sigma_r_avg(ie) = mean(diff_sr(:));
    R.err_sigma_t_max(ie) = max(diff_st(:));
    R.err_sigma_t_avg(ie) = mean(diff_st(:));

    fprintf('  ndof=%d  nel=%d  nlev=%d  wall=%.1fs\n', ...
            R.ndof(ie), R.nel(ie), R.nlevels(ie), R.wall(ie));
    fprintf('  err_u   max=%.4e  avg=%.4e\n', R.err_u_max(ie), R.err_u_avg(ie));
    fprintf('  sigma_r max=%.4e  avg=%.4e\n', R.err_sigma_r_max(ie), R.err_sigma_r_avg(ie));
    fprintf('  sigma_t max=%.4e  avg=%.4e\n', R.err_sigma_t_max(ie), R.err_sigma_t_avg(ie));
end

%% ========================================================================
% STEP 3: RESULTS SUMMARY
% =========================================================================
fprintf('\n\n');
fprintf('=====================================================================\n');
fprintf('  3D ROTATING CYLINDER: ESTIMATOR COMPARISON\n');
fprintf('=====================================================================\n');
fprintf('  Metric              |  div_sigma       |  stress_gradient\n');
fprintf('-------------------------------------------------------------\n');
fprintf('  DOFs                |  %10d      |  %10d\n', R.ndof(1), R.ndof(2));
fprintf('  Elements            |  %10d      |  %10d\n', R.nel(1), R.nel(2));
fprintf('  Levels              |  %10d      |  %10d\n', R.nlevels(1), R.nlevels(2));
fprintf('  Wall time [s]       |  %10.1f      |  %10.1f\n', R.wall(1), R.wall(2));
fprintf('  err_u max           |  %12.4e  |  %12.4e\n', R.err_u_max(1), R.err_u_max(2));
fprintf('  err_u avg           |  %12.4e  |  %12.4e\n', R.err_u_avg(1), R.err_u_avg(2));
fprintf('  err_sigma_r max     |  %12.4e  |  %12.4e\n', R.err_sigma_r_max(1), R.err_sigma_r_max(2));
fprintf('  err_sigma_r avg     |  %12.4e  |  %12.4e\n', R.err_sigma_r_avg(1), R.err_sigma_r_avg(2));
fprintf('  err_sigma_t max     |  %12.4e  |  %12.4e\n', R.err_sigma_t_max(1), R.err_sigma_t_max(2));
fprintf('  err_sigma_t avg     |  %12.4e  |  %12.4e\n', R.err_sigma_t_avg(1), R.err_sigma_t_avg(2));
fprintf('=====================================================================\n');

% Winner assessment
metrics = {'err_u_max','err_u_avg','err_sigma_r_max','err_sigma_r_avg','err_sigma_t_max','err_sigma_t_avg'};
wins_div = 0; wins_grad = 0;
for k = 1:numel(metrics)
    if R.(metrics{k})(1) < R.(metrics{k})(2)
        wins_div = wins_div + 1;
    else
        wins_grad = wins_grad + 1;
    end
end
fprintf('\n  div_sigma wins %d / %d metrics\n', wins_div, numel(metrics));
fprintf('  stress_gradient wins %d / %d metrics\n\n', wins_grad, numel(metrics));

% Comparison with 2D reference (plane-strain analytical)
fprintf('  Validation: 2D plane-strain analytical vs 3D numerical (last load step)\n');
fprintf('  (At full load, plastic response => analytical elastic is only a check\n');
fprintf('   for early steps; stress profiles will differ from elastic solution)\n');

%% ========================================================================
% STEP 4: PLOTS
% =========================================================================
radius_plot = linspace(R_i, R_o, n_eval);
est_colors = {[0.20 0.45 0.80], [0.85 0.40 0.20]};
est_styles = {'-o', '-s'};

% Load-displacement curves
figure(1); clf; hold on; grid on; box on;
plot(u_r_ref_vec, (0:nload)/nload, '-k', 'LineWidth', 1.5, 'DisplayName', 'Reference (3D)');
for ie = 1:numel(estimators)
    if ~isempty(R.u_r{ie})
        plot(R.u_r{ie}, (0:nload)/nload, est_styles{ie}, ...
             'Color', est_colors{ie}, 'LineWidth', 1.2, 'MarkerSize', 5, ...
             'MarkerFaceColor', est_colors{ie}, ...
             'DisplayName', strrep(estimators{ie}, '_', '\_'));
    end
end
xlabel('u_r at outer boundary [mm]');
ylabel('Load fraction');
title('3D Rotating cylinder: load-displacement');
legend('Location', 'northwest');

% Radial stress
figure(2); clf; hold on; grid on; box on;
for jload = 1:numel(load_steps_eval)
    plot(radius_plot, sigma_r_ref(:,jload), '-k', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Ref 3D, step %d', load_steps_eval(jload)));
end
for ie = 1:numel(estimators)
    if ~isempty(R.sigma_r{ie})
        for jload = 1:numel(load_steps_eval)
            plot(radius_plot, R.sigma_r{ie}(:,jload), est_styles{ie}, ...
                 'Color', est_colors{ie}, 'MarkerIndices', 1:5:n_eval, ...
                 'MarkerSize', 4, 'MarkerFaceColor', est_colors{ie}, ...
                 'DisplayName', sprintf('%s, step %d', ...
                     strrep(estimators{ie},'_','\_'), load_steps_eval(jload)));
        end
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_r [MPa]');
title('3D Rotating cylinder: radial stress');
legend('Location', 'best');

% Tangential stress
figure(3); clf; hold on; grid on; box on;
for jload = 1:numel(load_steps_eval)
    plot(radius_plot, sigma_t_ref(:,jload), '-k', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Ref 3D, step %d', load_steps_eval(jload)));
end
for ie = 1:numel(estimators)
    if ~isempty(R.sigma_t{ie})
        for jload = 1:numel(load_steps_eval)
            plot(radius_plot, R.sigma_t{ie}(:,jload), est_styles{ie}, ...
                 'Color', est_colors{ie}, 'MarkerIndices', 1:5:n_eval, ...
                 'MarkerSize', 4, 'MarkerFaceColor', est_colors{ie}, ...
                 'DisplayName', sprintf('%s, step %d', ...
                     strrep(estimators{ie},'_','\_'), load_steps_eval(jload)));
        end
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_t [MPa]');
title('3D Rotating cylinder: tangential stress');
legend('Location', 'best');

% Elastic validation
figure(4); clf; hold on; grid on; box on;
plot(r_ref, sig_r_el, '--r', 'LineWidth', 1.2, 'DisplayName', 'Analytical \sigma_r (plane strain)');
plot(r_ref, sig_t_el, '--b', 'LineWidth', 1.2, 'DisplayName', 'Analytical \sigma_t (plane strain)');
if ~isempty(R.sigma_r{1})
    plot(radius_plot, R.sigma_r{1}(:,1), '-ok', 'MarkerSize', 3, ...
         'MarkerIndices', 1:10:n_eval, 'DisplayName', sprintf('3D numerical, step %d', load_steps_eval(1)));
    plot(radius_plot, R.sigma_t{1}(:,1), '-sk', 'MarkerSize', 3, ...
         'MarkerIndices', 1:10:n_eval, 'DisplayName', sprintf('3D numerical \\sigma_t, step %d', load_steps_eval(1)));
end
xlabel('r [mm]'); ylabel('\sigma [MPa]');
title('3D: elastic validation against plane-strain analytical');
legend('Location', 'best');

saveas(figure(1), fullfile(results_dir, 'rot_cyl_3D_load_disp.png'));
saveas(figure(2), fullfile(results_dir, 'rot_cyl_3D_sigma_r.png'));
saveas(figure(3), fullfile(results_dir, 'rot_cyl_3D_sigma_t.png'));
saveas(figure(4), fullfile(results_dir, 'rot_cyl_3D_elastic_validation.png'));

%% ========================================================================
% STEP 5: SAVE
% =========================================================================
R.sigma_r_ref = sigma_r_ref;
R.sigma_t_ref = sigma_t_ref;
R.u_r_ref = u_r_ref_vec;
R.ref_ndof = cell_hspace_ref{end}.ndof;
R.ref_nel  = cell_hmsh_ref{end}.nel;
save(fullfile(results_dir, 'rotating_cylinder_3D_comparison.mat'), '-struct', 'R');

% CSV summary
fid = fopen(fullfile(results_dir, 'rotating_cylinder_3D_comparison.csv'), 'w');
fprintf(fid, 'estimator,ndof,nel,nlevels,wall_s,err_u_max,err_u_avg,err_sr_max,err_sr_avg,err_st_max,err_st_avg\n');
for ie = 1:numel(estimators)
    fprintf(fid, '%s,%d,%d,%d,%.3f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
            estimators{ie}, R.ndof(ie), R.nel(ie), R.nlevels(ie), R.wall(ie), ...
            R.err_u_max(ie), R.err_u_avg(ie), ...
            R.err_sigma_r_max(ie), R.err_sigma_r_avg(ie), ...
            R.err_sigma_t_max(ie), R.err_sigma_t_avg(ie));
end
fclose(fid);

fprintf('Results saved to %s\n', results_dir);
