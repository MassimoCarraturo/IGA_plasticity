% POSTPROCESS_QI_STUDY  Generate .dat tables and pgfplots from checkpoint data.
%
% Reads qi_study_checkpoint.mat and produces:
%   - convergence_<method>.dat
%   - load_disp_<method>.dat, load_disp_hill.dat
%   - sigma_r_step<N>.dat, sigma_t_step<N>.dat
%   - error_vs_load_<method>.dat
%   - qi_study_plots.tex  (standalone pgfplots document)
%
% Handles missing configurations (e.g. QI_p3 level 4 too expensive).

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
results_dir = fullfile(here, 'results', 'qi_study_v2');

addpath(genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'initialize')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));

cd(here);

%% Load checkpoint
tmp = load(fullfile(results_dir, 'qi_study_checkpoint.mat'), 'R');
R = tmp.R;

%% Study parameters (must match study script)
method_tags  = {'QI_C0', 'QI_p3', 'QI_graded'};
max_levels   = [1 2 3 4];
n_methods    = numel(method_tags);
n_levels     = numel(max_levels);
nload        = 5;
load_step_eval = [2 5];

% Hill analytical reference
E = 210000; nu = 0.3; sigma_y = 240;
a = 100; b = 200; P = 332 * 0.99;
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius_ex] = ...
    sphere_solution(E, nu, sigma_y, a, b, P, nload);
P_steps = linspace(P/nload, P, nload);

%% Determine which configs are completed
completed = false(n_methods, n_levels);
for im = 1:n_methods
    for il = 1:n_levels
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
        if isfield(R, tag) && R.(tag).ndof > 0
            completed(im, il) = true;
        end
    end
end

fprintf('Completed configurations:\n');
for im = 1:n_methods
    for il = 1:n_levels
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
        if completed(im, il)
            fprintf('  [OK]   %s  ndof=%d  err_t_rel=%.4e\n', tag, R.(tag).ndof, R.(tag).err_sigma_t_rel(end));
        else
            fprintf('  [MISS] %s\n', tag);
        end
    end
end

% Find best completed level per method (for stress profiles etc.)
best_level = zeros(n_methods, 1);
for im = 1:n_methods
    idx = find(completed(im,:), 1, 'last');
    if ~isempty(idx)
        best_level(im) = max_levels(idx);
    end
end
fprintf('\nBest level per method: ');
for im = 1:n_methods; fprintf('%s=%d  ', method_tags{im}, best_level(im)); end
fprintf('\n');

%% ========================================================================
%  CONSOLE SUMMARY TABLE
%  ========================================================================
fprintf('\n========================================================================\n');
fprintf('  QI VARIANT STUDY — SUMMARY (final load step)\n');
fprintf('========================================================================\n');
fprintf('  %-12s  %4s  %6s  %6s  %6s  %10s  %10s  %8s\n', ...
        'Config', 'Lev', 'DOFs', 'Nel', 'DOFs_s', 'err_sig_r', 'err_sig_t', 'Wall[s]');
fprintf('------------------------------------------------------------------------\n');
for im = 1:n_methods
    for il = 1:n_levels
        if ~completed(im,il); continue; end
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
        fprintf('  %-12s  %4d  %6d  %6d  %6d  %10.3e  %10.3e  %8.1f\n', ...
            method_tags{im}, max_levels(il), ...
            R.(tag).ndof, R.(tag).nel, R.(tag).ndof_scalar, ...
            R.(tag).err_sigma_r_rel(end), R.(tag).err_sigma_t_rel(end), ...
            R.(tag).wall);
    end
    fprintf('------------------------------------------------------------------------\n');
end

%% ========================================================================
%  EXPORT DATA TABLES
%  ========================================================================

% --- Table 1: convergence (only completed levels) -------------------------
for im = 1:n_methods
    fname = fullfile(results_dir, sprintf('convergence_%s.dat', method_tags{im}));
    fid = fopen(fname, 'w');
    fprintf(fid, 'max_level\tndof\tndof_scalar\tnel\terr_sigma_r\terr_sigma_t\terr_sigma_r_rel\terr_sigma_t_rel\twall\n');
    for il = 1:n_levels
        if ~completed(im,il); continue; end
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
        fprintf(fid, '%d\t%d\t%d\t%d\t%.10e\t%.10e\t%.10e\t%.10e\t%.4f\n', ...
            max_levels(il), R.(tag).ndof, R.(tag).ndof_scalar, R.(tag).nel, ...
            R.(tag).err_sigma_r(end), R.(tag).err_sigma_t(end), ...
            R.(tag).err_sigma_r_rel(end), R.(tag).err_sigma_t_rel(end), ...
            R.(tag).wall);
    end
    fclose(fid);
    fprintf('  Wrote %s\n', fname);
end

% --- Table 2: load-displacement (best level per method) -------------------
write_dat(fullfile(results_dir, 'load_disp_hill.dat'), {'u', 'P'}, [u_ex(:), P_ex(:)]);
for im = 1:n_methods
    tag = sprintf('%s_lev%d', method_tags{im}, best_level(im));
    write_dat(fullfile(results_dir, sprintf('load_disp_%s.dat', method_tags{im})), ...
              {'u', 'P'}, [R.(tag).u_r(:), R.(tag).P_i(:)]);
end

% --- Table 3: stress profiles at selected load steps ----------------------
for jl = 1:length(load_step_eval)
    step = load_step_eval(jl);

    % sigma_r
    headers = {'r', 'hill'};
    cols = [radius_ex(:), sigma_r_ex(step,:).'];
    for im = 1:n_methods
        tag = sprintf('%s_lev%d', method_tags{im}, best_level(im));
        headers{end+1} = method_tags{im}; %#ok<SAGROW>
        cols = [cols, R.(tag).sigma_rad(:, jl)]; %#ok<AGROW>
    end
    write_dat(fullfile(results_dir, sprintf('sigma_r_step%d.dat', step)), headers, cols);

    % sigma_t
    headers = {'r', 'hill'};
    cols = [radius_ex(:), sigma_t_ex(step,:).'];
    for im = 1:n_methods
        tag = sprintf('%s_lev%d', method_tags{im}, best_level(im));
        headers{end+1} = method_tags{im}; %#ok<SAGROW>
        cols = [cols, R.(tag).sigma_tan(:, jl)]; %#ok<AGROW>
    end
    write_dat(fullfile(results_dir, sprintf('sigma_t_step%d.dat', step)), headers, cols);
end

% --- Table 4: per-load-step error (best level per method) -----------------
for im = 1:n_methods
    tag = sprintf('%s_lev%d', method_tags{im}, best_level(im));
    fname = fullfile(results_dir, sprintf('error_vs_load_%s.dat', method_tags{im}));
    fid = fopen(fname, 'w');
    fprintf(fid, 'load_step\tP\terr_sigma_r\terr_sigma_t\terr_sigma_r_rel\terr_sigma_t_rel\n');
    for iload = 1:nload
        fprintf(fid, '%d\t%.6f\t%.10e\t%.10e\t%.10e\t%.10e\n', ...
            iload, P_steps(iload), ...
            R.(tag).err_sigma_r(iload), R.(tag).err_sigma_t(iload), ...
            R.(tag).err_sigma_r_rel(iload), R.(tag).err_sigma_t_rel(iload));
    end
    fclose(fid);
end

%% ========================================================================
%  GENERATE STANDALONE pgfplots DOCUMENT
%  ========================================================================
write_qi_pgfplots(fullfile(results_dir, 'qi_study_plots.tex'), ...
    method_tags, best_level, load_step_eval);

fprintf('\nAll done.  Results in %s\n', results_dir);

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function write_dat(filename, headers, data)
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', headers{1});
    for k = 2:numel(headers); fprintf(fid, '\t%s', headers{k}); end
    fprintf(fid, '\n');
    fmt = [strjoin(repmat({'%.10g'}, 1, size(data,2)), '\t') '\n'];
    for i = 1:size(data,1); fprintf(fid, fmt, data(i,:)); end
    fclose(fid);
    fprintf('  Wrote %s\n', filename);
end


function write_qi_pgfplots(filename, method_tags, best_level, load_steps)
    n_methods = numel(method_tags);

    % pgfplots style per method
    styles = {
        'green!50!black,  mark=diamond*,  mark size=2pt, very thick'
        'blue!70!black,   mark=*,         mark size=2pt, very thick'
        'red!70!black,    mark=triangle*, mark size=2pt, very thick'};
    labels = {'QI ($C^0$)', 'QI ($p{=}3$, $C^2$)', 'QI (graded)'};

    fid = fopen(filename, 'w');
    fp = @(varargin) fprintf(fid, varargin{:});

    fp('%% Auto-generated by postprocess_qi_study.m\n');
    fp('%% Compile: pdflatex qi_study_plots.tex\n');
    fp('\\documentclass[tikz,border=6pt]{standalone}\n');
    fp('\\usepackage{pgfplots}\n');
    fp('\\pgfplotsset{compat=1.18,\n');
    fp('  every axis/.append style={font=\\small, grid=both,\n');
    fp('    legend cell align=left}}\n');
    fp('\\usepgfplotslibrary{groupplots}\n\n');
    fp('\\begin{document}\n\n');

    % ==== Plot 1: Relative sigma_r error vs DOFs (log-log) ================
    fp('%%%% Relative sigma_r error vs DOFs\n');
    fp('\\begin{tikzpicture}\n');
    fp('\\begin{axis}[\n');
    fp('  width=11cm, height=8cm,\n');
    fp('  xlabel={Displacement DOFs},\n');
    fp('  ylabel={$\\|\\sigma_r - \\sigma_r^{\\mathrm{Hill}}\\|_{L^2} / \\|\\sigma_r^{\\mathrm{Hill}}\\|_{L^2}$},\n');
    fp('  title={Radial stress error vs.\\ refinement (final load step)},\n');
    fp('  xmode=log, ymode=log,\n');
    fp('  legend pos=north east]\n');
    for im = 1:n_methods
        fp('\\addplot[%s] table[x=ndof, y=err_sigma_r_rel] {convergence_%s.dat};\n', ...
           styles{im}, method_tags{im});
        fp('\\addlegendentry{%s}\n', labels{im});
    end
    fp('\\end{axis}\n');
    fp('\\end{tikzpicture}\n\n');

    % ==== Plot 2: Relative sigma_t error vs DOFs (log-log) ================
    fp('%%%% Relative sigma_t error vs DOFs\n');
    fp('\\begin{tikzpicture}\n');
    fp('\\begin{axis}[\n');
    fp('  width=11cm, height=8cm,\n');
    fp('  xlabel={Displacement DOFs},\n');
    fp('  ylabel={$\\|\\sigma_t - \\sigma_t^{\\mathrm{Hill}}\\|_{L^2} / \\|\\sigma_t^{\\mathrm{Hill}}\\|_{L^2}$},\n');
    fp('  title={Tangential stress error vs.\\ refinement (final load step)},\n');
    fp('  xmode=log, ymode=log,\n');
    fp('  legend pos=north east]\n');
    for im = 1:n_methods
        fp('\\addplot[%s] table[x=ndof, y=err_sigma_t_rel] {convergence_%s.dat};\n', ...
           styles{im}, method_tags{im});
        fp('\\addlegendentry{%s}\n', labels{im});
    end
    fp('\\end{axis}\n');
    fp('\\end{tikzpicture}\n\n');

    % ==== Plot 3: Load-displacement curves =================================
    fp('%%%% Load--displacement curves (best level per method)\n');
    fp('\\begin{tikzpicture}\n');
    fp('\\begin{axis}[\n');
    fp('  width=11cm, height=8cm,\n');
    fp('  xlabel={$u_r$ at outer boundary [mm]},\n');
    fp('  ylabel={Internal pressure $P$ [MPa]},\n');
    fp('  title={Load--displacement curve},\n');
    fp('  legend pos=north west]\n');
    fp('\\addplot[black, thick, no marks] table[x=u, y=P] {load_disp_hill.dat};\n');
    fp('\\addlegendentry{Hill (analytical)}\n');
    for im = 1:n_methods
        fp('\\addplot[%s] table[x=u, y=P] {load_disp_%s.dat};\n', ...
           styles{im}, method_tags{im});
        fp('\\addlegendentry{%s (lev.~%d)}\n', labels{im}, best_level(im));
    end
    fp('\\end{axis}\n');
    fp('\\end{tikzpicture}\n\n');

    % ==== Plots 4+5: Stress profiles at selected load steps ================
    for jl = 1:numel(load_steps)
        step = load_steps(jl);

        % sigma_r
        fp('%%%% Radial stress profile, load step %d\n', step);
        fp('\\begin{tikzpicture}\n');
        fp('\\begin{axis}[\n');
        fp('  width=11cm, height=8cm,\n');
        fp('  xlabel={Radial coordinate $r$ [mm]},\n');
        fp('  ylabel={$\\sigma_r$ [MPa]},\n');
        fp('  title={Radial stress, load step %d},\n', step);
        fp('  legend pos=south east]\n');
        fp('\\addplot[black, thick, no marks] table[x=r, y=hill] {sigma_r_step%d.dat};\n', step);
        fp('\\addlegendentry{Hill}\n');
        for im = 1:n_methods
            fp('\\addplot[%s, mark repeat=10] table[x=r, y=%s] {sigma_r_step%d.dat};\n', ...
               styles{im}, method_tags{im}, step);
            fp('\\addlegendentry{%s}\n', labels{im});
        end
        fp('\\end{axis}\n');
        fp('\\end{tikzpicture}\n\n');

        % sigma_t
        fp('%%%% Tangential stress profile, load step %d\n', step);
        fp('\\begin{tikzpicture}\n');
        fp('\\begin{axis}[\n');
        fp('  width=11cm, height=8cm,\n');
        fp('  xlabel={Radial coordinate $r$ [mm]},\n');
        fp('  ylabel={$\\sigma_t$ [MPa]},\n');
        fp('  title={Tangential stress, load step %d},\n', step);
        fp('  legend pos=south east]\n');
        fp('\\addplot[black, thick, no marks] table[x=r, y=hill] {sigma_t_step%d.dat};\n', step);
        fp('\\addlegendentry{Hill}\n');
        for im = 1:n_methods
            fp('\\addplot[%s, mark repeat=10] table[x=r, y=%s] {sigma_t_step%d.dat};\n', ...
               styles{im}, method_tags{im}, step);
            fp('\\addlegendentry{%s}\n', labels{im});
        end
        fp('\\end{axis}\n');
        fp('\\end{tikzpicture}\n\n');
    end

    % ==== Plot 6: Error vs load step (best level per method) ===============
    fp('%%%% Error evolution over load steps\n');
    fp('\\begin{tikzpicture}\n');
    fp('\\begin{groupplot}[\n');
    fp('  group style={group size=2 by 1, horizontal sep=2cm},\n');
    fp('  width=8cm, height=7cm,\n');
    fp('  xlabel={Load step},\n');
    fp('  ymode=log,\n');
    fp('  legend pos=north east]\n');

    % sigma_r vs load step
    fp('\\nextgroupplot[ylabel={$\\|\\sigma_r - \\sigma_r^{\\mathrm{Hill}}\\|_{L^2}$},\n');
    fp('  title={Radial stress}]\n');
    for im = 1:n_methods
        fp('\\addplot[%s] table[x=load_step, y=err_sigma_r] {error_vs_load_%s.dat};\n', ...
           styles{im}, method_tags{im});
        fp('\\addlegendentry{%s}\n', labels{im});
    end

    % sigma_t vs load step
    fp('\\nextgroupplot[ylabel={$\\|\\sigma_t - \\sigma_t^{\\mathrm{Hill}}\\|_{L^2}$},\n');
    fp('  title={Tangential stress}]\n');
    for im = 1:n_methods
        fp('\\addplot[%s] table[x=load_step, y=err_sigma_t] {error_vs_load_%s.dat};\n', ...
           styles{im}, method_tags{im});
    end

    fp('\\end{groupplot}\n');
    fp('\\end{tikzpicture}\n\n');

    fp('\\end{document}\n');
    fclose(fid);
    fprintf('  Wrote %s\n', filename);
end
