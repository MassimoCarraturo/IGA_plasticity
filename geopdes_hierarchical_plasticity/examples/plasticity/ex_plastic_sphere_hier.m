% EX_PLASTIC_SPHERE_HIER  Eighth-sphere J2-plasticity benchmark from
% "Computational Methods for Plasticity: Theory and Applications"
% (de Souza Neto, Peric, Owen).  Compares three projection methods for the
% history variables — L2, QI, QI_ref — against Hill's analytical solution.
%
% Outputs written next to this script:
%   sphere_results_<method>.mat     numerical curves for each method
%   pgfplots/sphere_load_disp.dat   load-displacement table (all methods + Hill)
%   pgfplots/sphere_sigma_r_<step>.dat   radial stress at chosen load step
%   pgfplots/sphere_sigma_t_<step>.dat   tangential stress at chosen load step
%   pgfplots/sphere_comparison.tex  standalone pgfplots document
%   sphere_comparison.fig / .png    MATLAB comparison figures

% 1) PHYSICAL DATA OF THE PROBLEM
clear all
clc
close all
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_eighth_sphere.txt'; % direction 1=radial

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [1];
problem_data.symm_sides   = [3 4 5];
problem_data.slider_sides = [6]; % zero displacement
problem_data.penalty_slider = @(x, y, z) 1e8 * ones (size (x));

% Physical parameters
E  =  210000;                                  % MPa
nu = 0.3;                                      % -
problem_data.yield_stress = @(x, y, z) 240 * ones (size (x));  % MPa
problem_data.kappa_lame = @(x, y, z) E/(3*(1-2*nu)) * ones (size (x));
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
nload = 5;
P = 332*.99;                                        % Limit internal pressure [MPa]
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.p = @(x, y, z) P*ones (size (x));

% Hill (analytical) reference
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution (E, nu, problem_data.yield_stress(1), 100, 200, P, nload);
load_step_eval_stress = [2, 5];

% 2) DISCRETIZATION COMMON TO ALL METHODS
p = 2;
method_data_template.degree     = [p p p];
method_data_template.regularity = [p-1 p-1 p-1];
method_data_template.nsub_coarse = [15, 5, 5];
method_data_template.nsub_refine = [2 2 2];
method_data_template.nquad      = [p+1 p+1 p+1];
method_data_template.space_type  = 'standard';
method_data_template.truncated   = 1;
method_data_template.nload      = nload;
method_data_template.newton_tol = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;

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

% 3) RUN EACH PROJECTION METHOD
proj_methods = {'L2', 'QI', 'QI_ref'};
results = struct();

for imethod = 1:numel(proj_methods)
    method_name = proj_methods{imethod};
    fprintf('\n================================================================\n');
    fprintf('Running ex_plastic_sphere_hier with type_projection = %s\n', method_name);
    fprintf('================================================================\n');

    method_data = method_data_template;
    method_data.type_projection = method_name;

    [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
     cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
        adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);

    % ---- Post-processing: load-displacement curve ---------------------------
    u_r = zeros(1, method_data.nload+1);
    P_i = zeros(1, method_data.nload+1);
    for i = 1:method_data.nload+1
        P_i(i) = P/method_data.nload * (i-1);
        if i > 1
            [eu, ~] = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1, .5, .5});
            u_r(i) = norm(eu);
        else
            u_r(i) = 0;
        end
    end

    % ---- Post-processing: radial / tangential stress profiles --------------
    rad_dir = [.5; .5; sqrt(2)/2];
    tan_dir = [.5; .5; -sqrt(2)/2];
    voigt   = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];
    rad_pos = linspace(0, 1, length(radius));

    sigma_rad = zeros(length(radius), length(load_step_eval_stress));
    sigma_tan = zeros(length(radius), length(load_step_eval_stress));

    for jload = 1:length(load_step_eval_stress)
        for jpoint = 1:length(rad_pos)
            position = {rad_pos(jpoint), 0.5, 0.75};
            sigma_matrix = zeros(3,3);
            for icomp = 1:6
                eu = sp_eval (cell_sigma{load_step_eval_stress(jload)}(:,icomp), ...
                              cell_hspace_scalar{load_step_eval_stress(jload)}, ...
                              geometry, position);
                sigma_matrix(voigt(icomp,1), voigt(icomp,2)) = eu;
                if icomp > 3
                    sigma_matrix(voigt(icomp,2), voigt(icomp,1)) = eu;
                end
            end
            sigma_rad(jpoint, jload) = rad_dir' * sigma_matrix * rad_dir;
            sigma_tan(jpoint, jload) = tan_dir' * sigma_matrix * tan_dir;
        end
    end

    R = struct();
    R.method      = method_name;
    R.u_r         = u_r;
    R.P_i         = P_i;
    R.radius      = radius;
    R.sigma_rad   = sigma_rad;
    R.sigma_tan   = sigma_tan;
    R.load_steps  = load_step_eval_stress;
    R.solution    = solution_data;
    results.(matlab.lang.makeValidName(method_name)) = R;

    save (sprintf('sphere_results_%s.mat', method_name), '-struct', 'R');
end

% 4) PLOT COMPARISON IN MATLAB
out_dir = 'pgfplots';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

method_styles = {'-o', '-s', '-d'};
method_colors = {[0.20 0.45 0.80], [0.85 0.40 0.20], [0.20 0.65 0.30]};
method_labels = proj_methods;

% Figure 1: load-displacement curve (Hill + three methods)
figure(1); clf; hold on; grid on; box on;
plot (u_ex, P_ex, '-k', 'LineWidth', 1.5, 'DisplayName', 'Hill (analytical)');
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    plot (R.u_r, R.P_i, method_styles{imethod}, ...
          'Color', method_colors{imethod}, 'LineWidth', 1.2, ...
          'MarkerSize', 6, 'MarkerFaceColor', method_colors{imethod}, ...
          'DisplayName', method_labels{imethod});
end
xlabel('u_r at outer boundary [mm]');
ylabel('Internal pressure P [MPa]');
title('Eighth-sphere: load–displacement curve');
legend('Location', 'northwest');

% Figure 2: radial stress profiles
figure(2); clf; hold on; grid on; box on;
for jload = 1:length(load_step_eval_stress)
    plot (radius, sigma_r_ex(load_step_eval_stress(jload),:), '-k', ...
          'LineWidth', 1.5, 'DisplayName', ...
          sprintf('Hill, step %d', load_step_eval_stress(jload)));
end
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    for jload = 1:length(load_step_eval_stress)
        plot (radius, R.sigma_rad(:, jload), method_styles{imethod}, ...
              'Color', method_colors{imethod}, ...
              'MarkerIndices', 1:10:length(radius), ...
              'DisplayName', sprintf('%s, step %d', ...
                                     method_labels{imethod}, ...
                                     load_step_eval_stress(jload)));
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_r [MPa]');
title('Eighth-sphere: radial stress');
legend('Location', 'best');

% Figure 3: tangential stress profiles
figure(3); clf; hold on; grid on; box on;
for jload = 1:length(load_step_eval_stress)
    plot (radius, sigma_t_ex(load_step_eval_stress(jload),:), '-k', ...
          'LineWidth', 1.5, 'DisplayName', ...
          sprintf('Hill, step %d', load_step_eval_stress(jload)));
end
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    for jload = 1:length(load_step_eval_stress)
        plot (radius, R.sigma_tan(:, jload), method_styles{imethod}, ...
              'Color', method_colors{imethod}, ...
              'MarkerIndices', 1:10:length(radius), ...
              'DisplayName', sprintf('%s, step %d', ...
                                     method_labels{imethod}, ...
                                     load_step_eval_stress(jload)));
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_t [MPa]');
title('Eighth-sphere: tangential stress');
legend('Location', 'best');

% Save MATLAB figures
saveas(figure(1), 'sphere_load_disp.png');
saveas(figure(2), 'sphere_sigma_r.png');
saveas(figure(3), 'sphere_sigma_t.png');

% 5) EXPORT DATA TABLES FOR pgfplots
% Load–displacement table (one column per method + Hill)
ld_table = cell(0,1);
ld_header = sprintf('u_hill\tP_hill');
ld_rows_hill = [u_ex(:), P_ex(:)];
% pad / merge: each method shares the same P_i (load steps known a priori),
% so we write a separate file per series for simplicity.
write_dat (fullfile(out_dir, 'sphere_load_disp_hill.dat'), ...
           {'u', 'P'}, [u_ex(:), P_ex(:)]);
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    write_dat (fullfile(out_dir, sprintf('sphere_load_disp_%s.dat', proj_methods{imethod})), ...
               {'u', 'P'}, [R.u_r(:), R.P_i(:)]);
end

% Stress profiles — one .dat per (quantity, load step), with a column per method + Hill
for jload = 1:length(load_step_eval_stress)
    step = load_step_eval_stress(jload);

    headers = {'r', 'hill'};
    cols    = [radius(:), sigma_r_ex(step,:).'];
    for imethod = 1:numel(proj_methods)
        R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
        headers{end+1} = proj_methods{imethod}; %#ok<SAGROW>
        cols = [cols, R.sigma_rad(:, jload)]; %#ok<AGROW>
    end
    write_dat (fullfile(out_dir, sprintf('sphere_sigma_r_step%d.dat', step)), headers, cols);

    headers = {'r', 'hill'};
    cols    = [radius(:), sigma_t_ex(step,:).'];
    for imethod = 1:numel(proj_methods)
        R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
        headers{end+1} = proj_methods{imethod}; %#ok<SAGROW>
        cols = [cols, R.sigma_tan(:, jload)]; %#ok<AGROW>
    end
    write_dat (fullfile(out_dir, sprintf('sphere_sigma_t_step%d.dat', step)), headers, cols);
end

% 6) GENERATE pgfplots LaTeX SCRIPT
write_pgfplots_tex (fullfile(out_dir, 'sphere_comparison.tex'), ...
                    proj_methods, load_step_eval_stress);

fprintf('\nAll done. Data tables and LaTeX in %s/\n', out_dir);

% =========================================================================
% Local helper functions
% =========================================================================
function write_dat (filename, headers, data)
% Write a tab-separated ASCII table for pgfplots `table[header=true]`.
    fid = fopen (filename, 'w');
    if fid < 0; error('Cannot open %s for writing', filename); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf (fid, '%s', headers{1});
    for k = 2:numel(headers); fprintf (fid, '\t%s', headers{k}); end
    fprintf (fid, '\n');
    ncols = size(data, 2);
    fmt = [strjoin(repmat({'%.10g'}, 1, ncols), '\t') '\n'];
    for i = 1:size(data,1)
        fprintf (fid, fmt, data(i,:));
    end
end

function write_pgfplots_tex (filename, proj_methods, load_steps)
    method_pgf_styles = {
        'mark=*,        mark size=2pt, color=blue!70!black,    very thick'
        'mark=square*,  mark size=2pt, color=orange!80!black,  very thick'
        'mark=diamond*, mark size=2pt, color=green!50!black,   very thick'};
    fid = fopen (filename, 'w');
    if fid < 0; error('Cannot open %s for writing', filename); end
    cleanup = onCleanup(@() fclose(fid));

    fprintf (fid, '%% Auto-generated by ex_plastic_sphere_hier.m\n');
    fprintf (fid, '%% Compile with: pdflatex sphere_comparison.tex\n');
    fprintf (fid, '\\documentclass[tikz,border=4pt]{standalone}\n');
    fprintf (fid, '\\usepackage{pgfplots}\n');
    fprintf (fid, '\\pgfplotsset{compat=1.18, every axis/.append style={font=\\small}}\n\n');
    fprintf (fid, '\\begin{document}\n\n');

    % --- Load-displacement plot ----------------------------------------------
    fprintf (fid, '\\begin{tikzpicture}\n');
    fprintf (fid, '\\begin{axis}[\n');
    fprintf (fid, '    width=10cm, height=7cm,\n');
    fprintf (fid, '    xlabel={$u_r$ at outer boundary [mm]},\n');
    fprintf (fid, '    ylabel={Internal pressure $P$ [MPa]},\n');
    fprintf (fid, '    title={Eighth-sphere: load--displacement curve},\n');
    fprintf (fid, '    legend pos=north west, grid=both]\n');
    fprintf (fid, '\\addplot[black, thick, no marks] table[x=u, y=P] {sphere_load_disp_hill.dat};\n');
    fprintf (fid, '\\addlegendentry{Hill (analytical)}\n');
    for k = 1:numel(proj_methods)
        fprintf (fid, '\\addplot[%s] table[x=u, y=P] {sphere_load_disp_%s.dat};\n', ...
                 method_pgf_styles{k}, proj_methods{k});
        fprintf (fid, '\\addlegendentry{%s}\n', latex_escape(proj_methods{k}));
    end
    fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');

    % --- Radial-stress plots --------------------------------------------------
    for jload = 1:numel(load_steps)
        step = load_steps(jload);
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{axis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={Radial coordinate $r$ [mm]},\n');
        fprintf (fid, '    ylabel={$\\sigma_r$ [MPa]},\n');
        fprintf (fid, '    title={Radial stress, load step %d},\n', step);
        fprintf (fid, '    legend pos=south east, grid=both]\n');
        fprintf (fid, '\\addplot[black, thick, no marks] table[x=r, y=hill] {sphere_sigma_r_step%d.dat};\n', step);
        fprintf (fid, '\\addlegendentry{Hill}\n');
        for k = 1:numel(proj_methods)
            fprintf (fid, '\\addplot[%s, mark repeat=10] table[x=r, y=%s] {sphere_sigma_r_step%d.dat};\n', ...
                     method_pgf_styles{k}, proj_methods{k}, step);
            fprintf (fid, '\\addlegendentry{%s}\n', latex_escape(proj_methods{k}));
        end
        fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    % --- Tangential-stress plots ---------------------------------------------
    for jload = 1:numel(load_steps)
        step = load_steps(jload);
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{axis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={Radial coordinate $r$ [mm]},\n');
        fprintf (fid, '    ylabel={$\\sigma_t$ [MPa]},\n');
        fprintf (fid, '    title={Tangential stress, load step %d},\n', step);
        fprintf (fid, '    legend pos=south east, grid=both]\n');
        fprintf (fid, '\\addplot[black, thick, no marks] table[x=r, y=hill] {sphere_sigma_t_step%d.dat};\n', step);
        fprintf (fid, '\\addlegendentry{Hill}\n');
        for k = 1:numel(proj_methods)
            fprintf (fid, '\\addplot[%s, mark repeat=10] table[x=r, y=%s] {sphere_sigma_t_step%d.dat};\n', ...
                     method_pgf_styles{k}, proj_methods{k}, step);
            fprintf (fid, '\\addlegendentry{%s}\n', latex_escape(proj_methods{k}));
        end
        fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    fprintf (fid, '\\end{document}\n');
end

function s = latex_escape(s)
    s = strrep(s, '_', '\_');
end
