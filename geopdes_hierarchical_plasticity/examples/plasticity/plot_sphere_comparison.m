% PLOT_SPHERE_COMPARISON  Load saved sphere_results_<method>.mat files,
% plot the L2/QI/QI_C0 comparison vs Hill's analytical solution, write
% pgfplots .dat tables, and emit a standalone pgfplots LaTeX script.
%
% Run this AFTER all three sphere_results_<method>.mat files exist
% (created by ex_plastic_sphere_hier or run_sphere_qiref_only).

here = fileparts(mfilename('fullpath'));
cd (here);
results_dir = fullfile(here, 'results');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

proj_methods = {'L2', 'QI', 'QI_C0'};

% Hill (analytical) reference — reuse the same constants the runs used
E  = 210000; nu = 0.3; sy = 240; nload = 5; P = 332*.99;
addpath (genpath(fullfile(here, '..', '..', '..', 'geopdes-3.2.2')));
addpath (genpath(fullfile(here, '..', '..', '..', 'nurbs-1.4.3')));
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution (E, nu, sy, 100, 200, P, nload);
load_step_eval_stress = [2, 5];

% Load all method results
results = struct();
for k = 1:numel(proj_methods)
    fname = fullfile(results_dir, sprintf('sphere_results_%s.mat', proj_methods{k}));
    if ~isfile(fname)
        error('Missing %s — run the corresponding solver first', fname);
    end
    R = load(fname);
    results.(matlab.lang.makeValidName(proj_methods{k})) = R;
end

% MATLAB plots
method_styles = {'-o', '-s', '-d'};
method_colors = {[0.20 0.45 0.80], [0.85 0.40 0.20], [0.20 0.65 0.30]};
method_labels = proj_methods;

figure(1); clf; hold on; grid on; box on;
plot (u_ex, P_ex, '-k', 'LineWidth', 1.5, 'DisplayName', 'Hill (analytical)');
for i = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{i}));
    plot (R.u_r, R.P_i, method_styles{i}, ...
          'Color', method_colors{i}, 'LineWidth', 1.2, ...
          'MarkerSize', 6, 'MarkerFaceColor', method_colors{i}, ...
          'DisplayName', method_labels{i});
end
xlabel('u_r at outer boundary [mm]');
ylabel('Internal pressure P [MPa]');
title('Eighth-sphere: load–displacement curve');
legend('Location', 'northwest');

figure(2); clf; hold on; grid on; box on;
for jload = 1:numel(load_step_eval_stress)
    plot (radius, sigma_r_ex(load_step_eval_stress(jload),:), '-k', ...
          'LineWidth', 1.5, 'DisplayName', ...
          sprintf('Hill, step %d', load_step_eval_stress(jload)));
end
for i = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{i}));
    for jload = 1:numel(load_step_eval_stress)
        plot (radius, R.sigma_rad(:, jload), method_styles{i}, ...
              'Color', method_colors{i}, ...
              'MarkerIndices', 1:10:length(radius), ...
              'DisplayName', sprintf('%s, step %d', method_labels{i}, ...
                                     load_step_eval_stress(jload)));
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_r [MPa]');
title('Eighth-sphere: radial stress');
legend('Location', 'best');

figure(3); clf; hold on; grid on; box on;
for jload = 1:numel(load_step_eval_stress)
    plot (radius, sigma_t_ex(load_step_eval_stress(jload),:), '-k', ...
          'LineWidth', 1.5, 'DisplayName', ...
          sprintf('Hill, step %d', load_step_eval_stress(jload)));
end
for i = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{i}));
    for jload = 1:numel(load_step_eval_stress)
        plot (radius, R.sigma_tan(:, jload), method_styles{i}, ...
              'Color', method_colors{i}, ...
              'MarkerIndices', 1:10:length(radius), ...
              'DisplayName', sprintf('%s, step %d', method_labels{i}, ...
                                     load_step_eval_stress(jload)));
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_t [MPa]');
title('Eighth-sphere: tangential stress');
legend('Location', 'best');

saveas(figure(1), fullfile(results_dir, 'sphere_load_disp.png'));
saveas(figure(2), fullfile(results_dir, 'sphere_sigma_r.png'));
saveas(figure(3), fullfile(results_dir, 'sphere_sigma_t.png'));

% pgfplots .dat tables
out_dir = fullfile(results_dir, 'pgfplots');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

write_dat (fullfile(out_dir, 'sphere_load_disp_hill.dat'), ...
           {'u', 'P'}, [u_ex(:), P_ex(:)]);
for i = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{i}));
    write_dat (fullfile(out_dir, sprintf('sphere_load_disp_%s.dat', proj_methods{i})), ...
               {'u', 'P'}, [R.u_r(:), R.P_i(:)]);
end

for jload = 1:numel(load_step_eval_stress)
    step = load_step_eval_stress(jload);

    headers = {'r', 'hill'};
    cols    = [radius(:), sigma_r_ex(step,:).'];
    for i = 1:numel(proj_methods)
        R = results.(matlab.lang.makeValidName(proj_methods{i}));
        headers{end+1} = proj_methods{i}; %#ok<SAGROW>
        cols = [cols, R.sigma_rad(:, jload)]; %#ok<AGROW>
    end
    write_dat (fullfile(out_dir, sprintf('sphere_sigma_r_step%d.dat', step)), headers, cols);

    headers = {'r', 'hill'};
    cols    = [radius(:), sigma_t_ex(step,:).'];
    for i = 1:numel(proj_methods)
        R = results.(matlab.lang.makeValidName(proj_methods{i}));
        headers{end+1} = proj_methods{i}; %#ok<SAGROW>
        cols = [cols, R.sigma_tan(:, jload)]; %#ok<AGROW>
    end
    write_dat (fullfile(out_dir, sprintf('sphere_sigma_t_step%d.dat', step)), headers, cols);
end

% pgfplots LaTeX
write_pgfplots_tex (fullfile(out_dir, 'sphere_comparison.tex'), ...
                    proj_methods, load_step_eval_stress);

fprintf('Plots saved as PNG; data tables and pgfplots LaTeX written to %s/\n', out_dir);

% ---------------- helpers ------------------------------------------------
function write_dat (filename, headers, data)
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

    fprintf (fid, '%% Auto-generated by plot_sphere_comparison.m\n');
    fprintf (fid, '%% Compile with:  pdflatex sphere_comparison.tex\n');
    fprintf (fid, '\\documentclass[tikz,border=4pt]{standalone}\n');
    fprintf (fid, '\\usepackage{pgfplots}\n');
    fprintf (fid, '\\pgfplotsset{compat=1.18, every axis/.append style={font=\\small}}\n\n');
    fprintf (fid, '\\begin{document}\n\n');

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
