% COMPARE_ESTIMATORS_ROTATING_CYLINDER  Demonstrate that the divergence-based
% (residual) error estimator outperforms the stress-gradient estimator when
% a non-trivial, spatially-varying body force is present, even in
% combination with a pressure-driven plastic front.
%
% Problem: thick-walled quarter-cylinder in plane strain loaded by:
%   (1) Internal pressure P at the inner radius — creates a sharp
%       elastic-plastic front that propagates outward (the classical
%       Hill/de Souza Neto benchmark).
%   (2) Centrifugal body force f = rho*omega^2 * [x; y] — adds a
%       spatially-varying equilibrium requirement across the domain.
%
% The combined loading is a realistic scenario (pressurised rotating
% turbine disk) where both a sharp plastic front AND a distributed body
% force coexist.
%
% WHY div_sigma wins here:
%   The stress-gradient estimator concentrates refinement at the sharp
%   elastic-plastic front (driven by pressure), but ignores the equilibrium
%   residual from the centrifugal body force in the rest of the domain.
%
%   The div_sigma estimator captures BOTH error sources:
%     - At the plastic front: div(sigma_h) changes rapidly → large residual
%     - In the elastic region: f + div(sigma_h) ≠ 0 where the coarse mesh
%       can't balance the spatially-varying body force
%   By distributing refinement optimally between both sources of error,
%   div_sigma achieves better overall accuracy per DOF.
%
% Reference solution: fine uniform mesh (no adaptivity).
%
% Output:
%   results/rotating_cylinder_comparison.mat
%   results/rotating_cylinder_comparison.csv

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

% --- Combined loading: internal pressure + centrifugal body force ---
%
% (1) Internal pressure: P drives a sharp plastic front from the inner
%     radius outward (Hill-type benchmark).
%     Elastic initial yield pressure for this geometry:
%       P_0 = sigma_y/sqrt(3) * (1 - R_i^2/R_o^2) = 103.9 MPa
%     We use P below P_0 — pressure alone does NOT cause yield.
P = 85;  % [MPa] internal pressure (~0.82 * P_0)

% (2) Centrifugal body force: f = rho*omega^2 * [x; y]
%     Choose rho*omega^2 so centrifugal stresses alone are ~60% of yield.
%     Neither loading alone causes plasticity, but combined they exceed
%     yield — a balanced scenario where both contributions matter equally.
rho_omega2 = 5.0e-3;  % [MPa/mm^2]

% Verify elastic stress levels (centrifugal only, for diagnostics)
[sig_r_el, sig_t_el, ~, r_ref] = rotating_cylinder_solution(E, nu, rho_omega2, R_i, R_o, 200);
sigma_eq_centrifugal = sqrt(3)/2 * abs(sig_t_el - sig_r_el);

% Elastic stress from pressure alone (at inner radius):
%   sigma_t_P = P*(R_o^2+R_i^2)/(R_o^2-R_i^2), sigma_r_P = -P
sigma_t_P = P * (R_o^2 + R_i^2) / (R_o^2 - R_i^2);
sigma_eq_pressure = sqrt(3)/2 * (sigma_t_P + P);  % von Mises at inner

fprintf('Combined loading check:\n');
fprintf('  Internal pressure P = %.1f MPa\n', P);
fprintf('  rho*omega^2 = %.4e MPa/mm^2\n', rho_omega2);
fprintf('  Centrifugal contribution at R_i: sigma_vM = %.1f MPa (%.0f%% of sigma_y)\n', ...
        sigma_eq_centrifugal(1), 100*sigma_eq_centrifugal(1)/sigma_y);
fprintf('  Pressure contribution at R_i:    sigma_vM = %.1f MPa (%.0f%% of sigma_y)\n', ...
        sigma_eq_pressure, 100*sigma_eq_pressure/sigma_y);
fprintf('  Combined (approx):               sigma_vM ~ %.1f MPa (%.0f%% of sigma_y)\n', ...
        sigma_eq_centrifugal(1) + sigma_eq_pressure, ...
        100*(sigma_eq_centrifugal(1) + sigma_eq_pressure)/sigma_y);
fprintf('  → Expect significant plasticity at inner radius\n');

%% ========================================================================
% PROBLEM DATA
% =========================================================================
nload = 10;  % load steps (ramp up pressure + angular velocity together)

problem_data.geo_name = 'geo_ring_SouzaNeto.txt';
problem_data.nmnn_sides   = [2];      % outer radius: traction-free (Neumann g=0)
problem_data.drchlt_sides = [];       % no Dirichlet
problem_data.press_sides  = [1];      % inner radius: pressure loading
problem_data.symm_sides   = [3 4];    % quarter symmetry
problem_data.slider_sides = [];

% Material
mu_lame    = E / (2*(1+nu));
kappa_lame = E / (3*(1-2*nu));
problem_data.yield_stress = @(x, y) sigma_y * ones(size(x));
problem_data.kappa_lame   = @(x, y) kappa_lame * ones(size(x));
problem_data.mu_lame      = @(x, y) mu_lame * ones(size(x));

% Body force: centrifugal f = rho*omega^2 * [x; y]
% The solver applies rhs_shape/nload * load_multiplier, so both pressure
% and body force are ramped linearly together over the load steps.
problem_data.f = @(x, y) rho_omega2 * cat(1, ...
    reshape(x, [1, size(x,1), size(x,2)]), ...
    reshape(y, [1, size(x,1), size(x,2)]));

% Pressure on inner surface (side 1)
problem_data.p = @(x, y) P * ones(size(x));

% Traction-free outer boundary
problem_data.g = @(x, y, ind) zeros(2, size(x,1), size(x,2));
problem_data.h = @(x, y, ind, mult) zeros(2, size(x,1), size(x,2));

%% ========================================================================
% DISCRETIZATION — common template
% =========================================================================
p = 2;
method_data_template.degree         = [p p];
method_data_template.regularity     = [p-1 p-1];
method_data_template.nsub_coarse    = [2 2];
method_data_template.nsub_refine    = [2 2];
method_data_template.nquad          = [p+1 p+1];
method_data_template.space_type     = 'standard';
method_data_template.truncated      = 1;
method_data_template.nload          = nload;
method_data_template.newton_tol     = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;
method_data_template.type_projection = 'QI_C0';  % will be overwritten

adaptivity_data_template.flag = 'elements';
adaptivity_data_template.C0_est = 1.0;
adaptivity_data_template.mark_param = 0.5;
adaptivity_data_template.mark_param_coarsening = 0.01;  % mild coarsening
adaptivity_data_template.mark_strategy = 'MS';
adaptivity_data_template.max_ndof  = 50000;
adaptivity_data_template.max_nel   = 20000;
adaptivity_data_template.num_max_iter = 3;   % 3 adaptive iterations per load step
adaptivity_data_template.tol = 1e-5 * 1e6;  % effectively never reached
adaptivity_data_template.adm_strategy = 'admissible';
adaptivity_data_template.coarsening_flag = 'any';
adaptivity_data_template.adm = p;
adaptivity_data_template.max_level = 5;

%% ========================================================================
% STEP 1: REFERENCE SOLUTION (fine uniform mesh, no adaptivity)
% =========================================================================
fprintf('\n================================================================\n');
fprintf('  REFERENCE SOLUTION (fine uniform mesh)\n');
fprintf('================================================================\n');

method_data_ref = method_data_template;
method_data_ref.nsub_coarse = [12 12];        % fine mesh for well-converged reference
method_data_ref.type_projection = 'QI';

adaptivity_data_ref = adaptivity_data_template;
adaptivity_data_ref.num_max_iter = 1;         % no adaptivity
adaptivity_data_ref.max_level = 1;
adaptivity_data_ref.mark_param_coarsening = 0; % disable coarsening for reference

wall0 = tic;
[geometry_ref, cell_hmsh_ref, cell_hspace_ref, cell_hspace_scalar_ref, ...
 cell_u_ref, ~, cell_sigma_ref, ~] = ...
    adaptivity_J2_plasticity(problem_data, method_data_ref, adaptivity_data_ref);
ref_time = toc(wall0);
fprintf('  Reference: ndof=%d, nel=%d, time=%.1fs\n', ...
        cell_hspace_ref{end}.ndof, cell_hmsh_ref{end}.nel, ref_time);

% Evaluate reference stress along radial line
n_eval = 200;
rad_pos = linspace(0, 1, n_eval);
pt_eval = {rad_pos, 0.5};  % parametric: xi1 = radial, xi2 = 0.5 (45° direction)

voigt_2d = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];  % for index mapping
load_steps_eval = [round(nload/2), nload];

sigma_r_ref = zeros(n_eval, numel(load_steps_eval));
sigma_t_ref = zeros(n_eval, numel(load_steps_eval));
u_r_ref_vec = zeros(1, nload+1);

% Reference displacement at outer boundary
for i = 1:nload
    [eu, ~] = sp_eval(cell_u_ref{i}, cell_hspace_ref{i}, geometry_ref, {1, 0.5});
    u_r_ref_vec(i+1) = norm(eu);
end

% Reference stress along radial line
for jload = 1:numel(load_steps_eval)
    ls = load_steps_eval(jload);
    for jpt = 1:n_eval
        position = {rad_pos(jpt), 0.5};
        sigma_xx = sp_eval(cell_sigma_ref{ls}(:,1), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        sigma_yy = sp_eval(cell_sigma_ref{ls}(:,2), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        sigma_xy = sp_eval(cell_sigma_ref{ls}(:,4), cell_hspace_scalar_ref{ls}, geometry_ref, position);
        % At theta = 45°: radial = [cos45, sin45], tangential = [-sin45, cos45]
        % sigma_rr = n^T * sigma * n
        c = sqrt(2)/2;
        sigma_r_ref(jpt, jload) = c^2*(sigma_xx + sigma_yy) + 2*c^2*sigma_xy;
        sigma_t_ref(jpt, jload) = c^2*(sigma_xx + sigma_yy) - 2*c^2*sigma_xy;
    end
end

%% ========================================================================
% STEP 2: ADAPTIVE SOLVES — div_sigma vs stress_gradient
% =========================================================================
estimators = {'div_sigma', 'stress_gradient'};
proj_method = 'QI';  % QI keeps hmsh_scalar == hmsh (no extra bisection)

R = struct();
R.estimators = estimators;
R.proj_method = proj_method;
R.rho_omega2 = rho_omega2;
R.load_steps_eval = load_steps_eval;
R.radius = linspace(R_i, R_o, n_eval);

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
    fprintf('  ADAPTIVE: estimator = %s,  method = %s\n', estimators{ie}, proj_method);
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
        [eu, ~] = sp_eval(cell_u{i}, cell_hspace{i}, geometry, {1, 0.5});
        u_r_ie(i+1) = norm(eu);
    end
    R.u_r{ie} = u_r_ie;

    % Evaluate stress along radial line
    sigma_r_ie = zeros(n_eval, numel(load_steps_eval));
    sigma_t_ie = zeros(n_eval, numel(load_steps_eval));
    for jload = 1:numel(load_steps_eval)
        ls = load_steps_eval(jload);
        for jpt = 1:n_eval
            position = {rad_pos(jpt), 0.5};
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
fprintf('╔═══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  PRESSURISED ROTATING CYLINDER: ESTIMATOR COMPARISON                 ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Metric              │  div_sigma       │  stress_gradient           ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  DOFs                │  %10d      │  %10d               ║\n', R.ndof(1), R.ndof(2));
fprintf('║  Elements            │  %10d      │  %10d               ║\n', R.nel(1), R.nel(2));
fprintf('║  Wall time [s]       │  %10.1f      │  %10.1f               ║\n', R.wall(1), R.wall(2));
fprintf('║  err_u max           │  %12.4e  │  %12.4e           ║\n', R.err_u_max(1), R.err_u_max(2));
fprintf('║  err_u avg           │  %12.4e  │  %12.4e           ║\n', R.err_u_avg(1), R.err_u_avg(2));
fprintf('║  err_sigma_r max     │  %12.4e  │  %12.4e           ║\n', R.err_sigma_r_max(1), R.err_sigma_r_max(2));
fprintf('║  err_sigma_r avg     │  %12.4e  │  %12.4e           ║\n', R.err_sigma_r_avg(1), R.err_sigma_r_avg(2));
fprintf('║  err_sigma_t max     │  %12.4e  │  %12.4e           ║\n', R.err_sigma_t_max(1), R.err_sigma_t_max(2));
fprintf('║  err_sigma_t avg     │  %12.4e  │  %12.4e           ║\n', R.err_sigma_t_avg(1), R.err_sigma_t_avg(2));
fprintf('╚═══════════════════════════════════════════════════════════════════════╝\n');

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

%% ========================================================================
% STEP 4: PLOTS
% =========================================================================
radius_plot = linspace(R_i, R_o, n_eval);
est_colors = {[0.20 0.45 0.80], [0.85 0.40 0.20]};
est_styles = {'-o', '-s'};

% Load-displacement curves
figure(1); clf; hold on; grid on; box on;
plot(u_r_ref_vec, (0:nload)/nload, '-k', 'LineWidth', 1.5, 'DisplayName', 'Reference');
for ie = 1:numel(estimators)
    if ~isempty(R.u_r{ie})
        plot(R.u_r{ie}, (0:nload)/nload, est_styles{ie}, ...
             'Color', est_colors{ie}, 'LineWidth', 1.2, 'MarkerSize', 5, ...
             'MarkerFaceColor', est_colors{ie}, ...
             'DisplayName', strrep(estimators{ie}, '_', '\_'));
    end
end
xlabel('u_r at outer boundary [mm]');
ylabel('Load fraction (P + \omega^2) / max');
title('Pressurised rotating cylinder: load-displacement');
legend('Location', 'northwest');

% Radial stress
figure(2); clf; hold on; grid on; box on;
for jload = 1:numel(load_steps_eval)
    plot(radius_plot, sigma_r_ref(:,jload), '-k', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Ref, step %d', load_steps_eval(jload)));
end
for ie = 1:numel(estimators)
    if ~isempty(R.sigma_r{ie})
        for jload = 1:numel(load_steps_eval)
            plot(radius_plot, R.sigma_r{ie}(:,jload), est_styles{ie}, ...
                 'Color', est_colors{ie}, 'MarkerIndices', 1:10:n_eval, ...
                 'MarkerSize', 4, 'MarkerFaceColor', est_colors{ie}, ...
                 'DisplayName', sprintf('%s, step %d', ...
                     strrep(estimators{ie},'_','\_'), load_steps_eval(jload)));
        end
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_r [MPa]');
title('Pressurised rotating cylinder: radial stress');
legend('Location', 'best');

% Tangential stress
figure(3); clf; hold on; grid on; box on;
for jload = 1:numel(load_steps_eval)
    plot(radius_plot, sigma_t_ref(:,jload), '-k', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Ref, step %d', load_steps_eval(jload)));
end
for ie = 1:numel(estimators)
    if ~isempty(R.sigma_t{ie})
        for jload = 1:numel(load_steps_eval)
            plot(radius_plot, R.sigma_t{ie}(:,jload), est_styles{ie}, ...
                 'Color', est_colors{ie}, 'MarkerIndices', 1:10:n_eval, ...
                 'MarkerSize', 4, 'MarkerFaceColor', est_colors{ie}, ...
                 'DisplayName', sprintf('%s, step %d', ...
                     strrep(estimators{ie},'_','\_'), load_steps_eval(jload)));
        end
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_t [MPa]');
title('Pressurised rotating cylinder: tangential stress');
legend('Location', 'best');

% Elastic analytical vs reference (validation)
figure(4); clf; hold on; grid on; box on;
plot(r_ref, sig_r_el, '--r', 'LineWidth', 1.2, 'DisplayName', 'Elastic analytical \sigma_r');
plot(r_ref, sig_t_el, '--b', 'LineWidth', 1.2, 'DisplayName', 'Elastic analytical \sigma_t');
% Plot first load step (should be elastic)
if ~isempty(R.sigma_r{1})
    plot(radius_plot, R.sigma_r{1}(:,1), '-ok', 'MarkerSize', 3, ...
         'MarkerIndices', 1:20:n_eval, 'DisplayName', sprintf('Numerical, step %d', load_steps_eval(1)));
end
xlabel('r [mm]'); ylabel('\sigma [MPa]');
title('Elastic validation: analytical vs numerical');
legend('Location', 'best');

saveas(figure(1), fullfile(results_dir, 'rot_cyl_load_disp.png'));
saveas(figure(2), fullfile(results_dir, 'rot_cyl_sigma_r.png'));
saveas(figure(3), fullfile(results_dir, 'rot_cyl_sigma_t.png'));
saveas(figure(4), fullfile(results_dir, 'rot_cyl_elastic_validation.png'));

%% ========================================================================
% STEP 5: SAVE
% =========================================================================
R.sigma_r_ref = sigma_r_ref;
R.sigma_t_ref = sigma_t_ref;
R.u_r_ref = u_r_ref_vec;
R.ref_ndof = cell_hspace_ref{end}.ndof;
R.ref_nel  = cell_hmsh_ref{end}.nel;
save(fullfile(results_dir, 'rotating_cylinder_comparison.mat'), '-struct', 'R');

% CSV summary
fid = fopen(fullfile(results_dir, 'rotating_cylinder_comparison.csv'), 'w');
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
