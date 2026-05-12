% ADAPTIVITY_MAX_LEVEL_STUDY  Sweep adaptivity_data.max_level from 1 to 4
% on the eighth-sphere J2-plasticity benchmark, starting from a coarser
% initial mesh nsub_coarse = [3 1 1] and allowing num_max_iter = 5
% refinement loops per load step.  Repeats the sweep for the three
% projection strategies L2 / QI / QI_C0 and reports both the average and
% maximum error against Hill's analytical solution at the chosen load
% steps.  Writes a standalone pgfplots LaTeX document plus the underlying
% .dat tables.
%
% Supports resume: if adaptivity_study.mat already exists, configs whose
% results are present (non-NaN wallclock) are skipped.
%
% Outputs (in this folder):
%   adaptivity_study.mat           per-(method, level) results struct
%   adaptivity_study_summary.csv   one row per (method, level)
%   pgfplots/adaptivity_summary.dat            master table for pgfplots
%   pgfplots/adaptivity_load_disp_<m>_lvl<L>.dat  per-method, per-level load-disp curves
%   pgfplots/adaptivity_sigma_r_<m>_lvl<L>.dat    per-method, per-level radial stress
%   pgfplots/adaptivity_sigma_t_<m>_lvl<L>.dat    per-method, per-level tangential stress
%   pgfplots/adaptivity_study.tex              standalone pgfplots document

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');

% ---- MATLAB path setup --------------------------------------------------
addpath (genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath (genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'initialize')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
addpath (genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
addpath (fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

cd (here);
warning ('off', 'MATLAB:nearlySingularMatrix');
warning ('off', 'MATLAB:singularMatrix');

% ---- Problem data (same physics as ex_plastic_sphere_hier) --------------
problem_data.geo_name = 'geo_eighth_sphere.txt';
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [1];
problem_data.symm_sides   = [3 4 5];
problem_data.slider_sides = [6];
problem_data.penalty_slider = @(x, y, z) 1e8 * ones (size (x));

E  = 210000; nu = 0.3;
problem_data.yield_stress = @(x, y, z) 240 * ones (size (x));
problem_data.kappa_lame   = @(x, y, z) E/(3*(1-2*nu)) * ones (size (x));
problem_data.mu_lame      = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

nload = 5;
P     = 332*.99;
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.p = @(x, y, z) P*ones (size (x));

% Hill (analytical) reference
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution (E, nu, problem_data.yield_stress(1), 100, 200, P, nload);
load_step_eval_stress = [2, 5];

% ---- Discretization (coarser initial mesh, more adaptive iterations) ----
p = 2;
method_data_template.degree         = [p p p];
method_data_template.regularity     = [p-1 p-1 p-1];
method_data_template.nsub_coarse    = [3, 1, 1];
method_data_template.nsub_refine    = [2 2 2];
method_data_template.nquad          = [p+1 p+1 p+1];
method_data_template.space_type     = 'standard';
method_data_template.truncated      = 1;
method_data_template.nload          = nload;
method_data_template.newton_tol     = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;

adaptivity_data_template.flag             = 'elements';
adaptivity_data_template.C0_est           = 1.0;
adaptivity_data_template.mark_param       = 0.9;
adaptivity_data_template.mark_param_coarsening = 0.05;
adaptivity_data_template.mark_strategy    = 'MS';
adaptivity_data_template.max_ndof         = 200000;
adaptivity_data_template.max_nel          = 100000;
adaptivity_data_template.num_max_iter     = 5;
adaptivity_data_template.tol              = 1e-5 * 3.14 * 30000;
adaptivity_data_template.adm_strategy     = 'admissible';
adaptivity_data_template.coarsening_flag  = 'any';
adaptivity_data_template.adm              = p;

% ---- Sweep over (method, max_level) -------------------------------------
methods    = {'L2', 'QI', 'QI_C0'};
max_levels = 1:4;
nM = numel(methods);
nL = numel(max_levels);

% ---- Resume: load existing results if available -------------------------
mat_file = fullfile(here, 'adaptivity_study.mat');
if exist(mat_file, 'file')
    R_old = load(mat_file);
    nL_old = numel(R_old.max_levels);
    nM_old = numel(R_old.methods);
    R = struct();
    R.methods    = methods;
    R.max_levels = max_levels;
    R.radius     = radius;
    R.load_steps = load_step_eval_stress;
    for fld = {'wallclock_s','cpu_s','ndof_final','nel_final','nlevels_final', ...
               'err_u_inf','err_u_avg','err_u_rel', ...
               'err_sigma_r_inf','err_sigma_r_avg', ...
               'err_sigma_t_inf','err_sigma_t_avg'}
        R.(fld{1}) = nan(nM, nL);
        if isfield(R_old, fld{1})
            R.(fld{1})(1:nM_old, 1:nL_old) = R_old.(fld{1});
        end
    end
    R.u_r_curves       = cell(nM, nL);
    R.P_i_curves       = cell(nM, nL);
    R.sigma_rad_curves = cell(nM, nL);
    R.sigma_tan_curves = cell(nM, nL);
    for fld = {'u_r_curves','P_i_curves','sigma_rad_curves','sigma_tan_curves'}
        if isfield(R_old, fld{1})
            R.(fld{1})(1:nM_old, 1:nL_old) = R_old.(fld{1});
        end
    end
    fprintf('Loaded existing results (%d methods x %d levels)\n', nM_old, nL_old);
else
    R = struct();
    R.methods    = methods;
    R.max_levels = max_levels;
    R.radius     = radius;
    R.load_steps = load_step_eval_stress;
    for fld = {'wallclock_s','cpu_s','ndof_final','nel_final','nlevels_final', ...
               'err_u_inf','err_u_avg','err_u_rel', ...
               'err_sigma_r_inf','err_sigma_r_avg', ...
               'err_sigma_t_inf','err_sigma_t_avg'}
        R.(fld{1}) = nan(nM, nL);
    end
    R.u_r_curves      = cell(nM, nL);
    R.P_i_curves      = cell(nM, nL);
    R.sigma_rad_curves = cell(nM, nL);
    R.sigma_tan_curves = cell(nM, nL);
end

for im = 1:nM
    for il = 1:nL
        method_name = methods{im};
        ml = max_levels(il);

        % Skip already-computed configs (resume support)
        if ~isnan(R.wallclock_s(im, il))
            fprintf('  %-7s  max_level = %d   SKIPPED (already computed)\n', method_name, ml);
            continue;
        end

        fprintf('\n========================================================\n');
        fprintf('  %-7s  max_level = %d   (sweep %d / %d)\n', method_name, ml, ...
                (im-1)*nL + il, nM*nL);
        fprintf('========================================================\n');

        method_data     = method_data_template;
        method_data.type_projection = method_name;
        adaptivity_data = adaptivity_data_template;
        adaptivity_data.max_level = ml;

        cpu0  = cputime();
        wall0 = tic;
        try
            [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
             cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
                adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);
        catch e
            fprintf('  FAILED: %s\n', e.message);
            R.wallclock_s(im, il) = toc(wall0);
            R.cpu_s(im, il)       = cputime() - cpu0;
            save (fullfile(here, 'adaptivity_study.mat'), '-struct', 'R');
            continue;
        end
        R.wallclock_s(im, il) = toc(wall0);
        R.cpu_s(im, il)       = cputime() - cpu0;

        final_hspace = cell_hspace{end};
        final_hmsh   = cell_hmsh{end};
        R.ndof_final(im, il)    = final_hspace.ndof;
        R.nel_final(im, il)     = final_hmsh.nel;
        R.nlevels_final(im, il) = final_hmsh.nlevels;

        % --- load-displacement curve --------------------------------------------
        u_r = zeros(1, nload+1);
        P_i = zeros(1, nload+1);
        for i = 1:nload+1
            P_i(i) = P/nload * (i-1);
            if i > 1
                [eu, ~] = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1, .5, .5});
                u_r(i) = norm(eu);
            end
        end
        R.u_r_curves{im, il} = u_r;
        R.P_i_curves{im, il} = P_i;

        Pi_pos    = P_i(P_i > 0);
        u_r_pos   = u_r(P_i > 0);
        keep      = [true, diff(P_ex) > 0];
        hill_u_at = interp1(P_ex(keep), u_ex(keep), Pi_pos, 'linear', 'extrap');
        diff_u    = abs(u_r_pos - hill_u_at);
        R.err_u_inf(im, il) = max(diff_u);
        R.err_u_avg(im, il) = mean(diff_u);
        R.err_u_rel(im, il) = R.err_u_inf(im, il) / max(abs(hill_u_at));

        % --- stress profiles at chosen load steps --------------------------------
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
        R.sigma_rad_curves{im, il} = sigma_rad;
        R.sigma_tan_curves{im, il} = sigma_tan;

        % Errors against Hill at all evaluated load steps (max & avg)
        diff_sr = zeros(length(radius), length(load_step_eval_stress));
        diff_st = zeros(length(radius), length(load_step_eval_stress));
        for jload = 1:length(load_step_eval_stress)
            ls = load_step_eval_stress(jload);
            diff_sr(:, jload) = abs(sigma_rad(:, jload) - sigma_r_ex(ls, :).');
            diff_st(:, jload) = abs(sigma_tan(:, jload) - sigma_t_ex(ls, :).');
        end
        R.err_sigma_r_inf(im, il) = max(diff_sr(:));
        R.err_sigma_r_avg(im, il) = mean(diff_sr(:));
        R.err_sigma_t_inf(im, il) = max(diff_st(:));
        R.err_sigma_t_avg(im, il) = mean(diff_st(:));

        fprintf('  ndof_final = %d   nel_final = %d   nlevels = %d\n', ...
                R.ndof_final(im, il), R.nel_final(im, il), R.nlevels_final(im, il));
        fprintf('  wall = %.2f s   CPU = %.2f s\n', R.wallclock_s(im, il), R.cpu_s(im, il));
        fprintf('  err_u   max=%.3e  avg=%.3e   sigma_r max=%.3e avg=%.3e   sigma_t max=%.3e avg=%.3e\n', ...
                R.err_u_inf(im, il), R.err_u_avg(im, il), ...
                R.err_sigma_r_inf(im, il), R.err_sigma_r_avg(im, il), ...
                R.err_sigma_t_inf(im, il), R.err_sigma_t_avg(im, il));

        save (fullfile(here, 'adaptivity_study.mat'), '-struct', 'R');
    end
end

% ---- CSV summary --------------------------------------------------------
fid = fopen (fullfile(here, 'adaptivity_study_summary.csv'), 'w');
cleanup_csv = onCleanup(@() fclose(fid));
fprintf (fid, 'method,max_level,ndof_final,nel_final,nlevels_final,wallclock_s,cpu_s,err_u_max,err_u_avg,err_sigma_r_max,err_sigma_r_avg,err_sigma_t_max,err_sigma_t_avg\n');
for im = 1:nM
    for il = 1:nL
        fprintf (fid, '%s,%d,%d,%d,%d,%.3f,%.3f,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n', ...
                 methods{im}, max_levels(il), ...
                 R.ndof_final(im,il), R.nel_final(im,il), R.nlevels_final(im,il), ...
                 R.wallclock_s(im,il), R.cpu_s(im,il), ...
                 R.err_u_inf(im,il), R.err_u_avg(im,il), ...
                 R.err_sigma_r_inf(im,il), R.err_sigma_r_avg(im,il), ...
                 R.err_sigma_t_inf(im,il), R.err_sigma_t_avg(im,il));
    end
end
clear cleanup_csv;
fprintf('Wrote adaptivity_study.mat and adaptivity_study_summary.csv\n');

% ---- pgfplots .dat tables -----------------------------------------------
out_dir = 'pgfplots';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% Master summary, one .dat per method, with one row per max_level
for im = 1:nM
    rows = [max_levels(:), R.ndof_final(im,:).', R.nel_final(im,:).', ...
            R.wallclock_s(im,:).', R.cpu_s(im,:).', ...
            R.err_u_inf(im,:).', R.err_u_avg(im,:).', ...
            R.err_sigma_r_inf(im,:).', R.err_sigma_r_avg(im,:).', ...
            R.err_sigma_t_inf(im,:).', R.err_sigma_t_avg(im,:).'];
    write_dat (fullfile(out_dir, sprintf('adaptivity_summary_%s.dat', methods{im})), ...
               {'max_level','ndof','nel','wallclock','cpu', ...
                'err_u_max','err_u_avg','err_sigma_r_max','err_sigma_r_avg', ...
                'err_sigma_t_max','err_sigma_t_avg'}, rows);
end

% Hill reference curves
write_dat (fullfile(out_dir, 'adaptivity_load_disp_hill.dat'), ...
           {'u','P'}, [u_ex(:), P_ex(:)]);
last_step = load_step_eval_stress(end);
write_dat (fullfile(out_dir, sprintf('adaptivity_sigma_r_hill_step%d.dat', last_step)), ...
           {'r','hill'}, [radius(:), sigma_r_ex(last_step,:).']);
write_dat (fullfile(out_dir, sprintf('adaptivity_sigma_t_hill_step%d.dat', last_step)), ...
           {'r','hill'}, [radius(:), sigma_t_ex(last_step,:).']);

% Per-(method, level) curves
for im = 1:nM
    for il = 1:nL
        if isempty(R.u_r_curves{im, il}); continue; end
        ml = max_levels(il);
        m  = methods{im};
        write_dat (fullfile(out_dir, sprintf('adaptivity_load_disp_%s_lvl%d.dat', m, ml)), ...
                   {'u','P'}, [R.u_r_curves{im, il}(:), R.P_i_curves{im, il}(:)]);
        write_dat (fullfile(out_dir, sprintf('adaptivity_sigma_r_%s_lvl%d.dat', m, ml)), ...
                   {'r','sigma'}, [radius(:), R.sigma_rad_curves{im, il}(:, end)]);
        write_dat (fullfile(out_dir, sprintf('adaptivity_sigma_t_%s_lvl%d.dat', m, ml)), ...
                   {'r','sigma'}, [radius(:), R.sigma_tan_curves{im, il}(:, end)]);
    end
end

% ---- pgfplots LaTeX document --------------------------------------------
write_pgfplots_tex (fullfile(out_dir, 'adaptivity_study.tex'), ...
                    methods, max_levels, last_step);

fprintf('\nAll outputs written.  pgfplots master file: %s\n', ...
        fullfile(out_dir, 'adaptivity_study.tex'));

% ---- MATLAB plots --------------------------------------------------------
plot_adaptivity_results (R, methods, max_levels, u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius, last_step);

% =========================================================================
% Local helpers
% =========================================================================
function write_dat (filename, headers, data)
    fid = fopen (filename, 'w');
    if fid < 0; error('Cannot open %s', filename); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf (fid, '%s', headers{1});
    for k = 2:numel(headers); fprintf (fid, '\t%s', headers{k}); end
    fprintf (fid, '\n');
    % Build the format string explicitly: single-quoted '\t' / '\n' are
    % 2-char literals, so the previous `fmt(1:end-1)` truncation left a
    % stray backslash and turned the row separator into the literal "\n".
    fmt = [repmat('%.10g\t', 1, size(data, 2) - 1), '%.10g\n'];
    for i = 1:size(data,1)
        fprintf (fid, fmt, data(i,:));
    end
end

function write_pgfplots_tex (filename, methods, max_levels, last_step)
    method_colors = {'blue!70!black', 'orange!80!black', 'green!50!black'};
    method_marks  = {'*', 'square*', 'triangle*'};
    % Legend labels matching user conventions: L2 -> $L_2$, QI -> QI $C^1$, QI_C0 -> QI $C^0$
    method_legend = containers.Map( ...
        {'L2', 'QI', 'QI_C0'}, ...
        {'$L_2$', 'QI $C^1$', 'QI $C^0$'});
    level_colors  = {'blue!70!black', 'orange!80!black', 'green!50!black', 'red!70!black'};
    level_marks   = {'*', 'square*', 'triangle*', 'diamond*'};

    xtick_str = ['{', strjoin(arrayfun(@num2str, max_levels, 'UniformOutput', false), ','), '}'];

    fid = fopen (filename, 'w');
    if fid < 0; error('Cannot open %s', filename); end
    cleanup = onCleanup(@() fclose(fid));

    fprintf (fid, '%% Auto-generated by adaptivity_max_level_study.m\n');
    fprintf (fid, '%% Compile with:\n');
    fprintf (fid, '%%   pdflatex adaptivity_study.tex\n');
    fprintf (fid, '%%   magick -density 300 adaptivity_study.pdf adaptivity_study.png\n');
    fprintf (fid, '%%   (produces adaptivity_study-0.png .. adaptivity_study-N.png)\n');
    fprintf (fid, '\\documentclass[tikz,border=4pt,convert=true,magick]{standalone}\n');
    fprintf (fid, '\\usepackage{pgfplots}\n');
    fprintf (fid, '\\pgfplotsset{compat=1.18, every axis/.append style={font=\\small}}\n\n');
    fprintf (fid, '\\begin{document}\n\n');

    % ---- Convergence plots (one per error metric, lines = methods) -----
    metrics = {
        'err_u_max',       '$\Vert u_r - u_\mathrm{Hill}\Vert_\infty$',  'Max u-error';
        'err_u_avg',       '$\overline{|u_r - u_\mathrm{Hill}|}$',        'Avg u-error';
        'err_sigma_r_max', '$\Vert\sigma_r - \sigma_r^\mathrm{Hill}\Vert_\infty$', 'Max $\sigma_r$ error';
        'err_sigma_r_avg', '$\overline{|\sigma_r - \sigma_r^\mathrm{Hill}|}$',     'Avg $\sigma_r$ error';
        'err_sigma_t_max', '$\Vert\sigma_t - \sigma_t^\mathrm{Hill}\Vert_\infty$', 'Max $\sigma_t$ error';
        'err_sigma_t_avg', '$\overline{|\sigma_t - \sigma_t^\mathrm{Hill}|}$',     'Avg $\sigma_t$ error';
    };
    for k = 1:size(metrics,1)
        col   = metrics{k,1};
        ylab  = metrics{k,2};
        title = metrics{k,3};
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{semilogyaxis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={Maximum number of levels},\n');
        fprintf (fid, '    ylabel={%s},\n', ylab);
        fprintf (fid, '    title={%s vs max\\_level},\n', title);
        fprintf (fid, '    legend pos=outer north east, grid=both, xtick=%s]\n', xtick_str);
        for im = 1:numel(methods)
            fprintf (fid, '\\addplot[%s, mark=%s, mark size=2.5pt, very thick] table[x=max_level, y=%s] {adaptivity_summary_%s.dat};\n', ...
                     method_colors{1+mod(im-1,numel(method_colors))}, ...
                     method_marks {1+mod(im-1,numel(method_marks))},  ...
                     col, methods{im});
            fprintf (fid, '\\addlegendentry{%s}\n', method_legend(methods{im}));
        end
        fprintf (fid, '\\end{semilogyaxis}\n\\end{tikzpicture}\n\n');
    end

    % ---- Computational cost vs max_level (per method) -------------------
    fprintf (fid, '\\begin{tikzpicture}\n');
    fprintf (fid, '\\begin{axis}[\n');
    fprintf (fid, '    width=10cm, height=7cm,\n');
    fprintf (fid, '    xlabel={Maximum number of levels},\n');
    fprintf (fid, '    ylabel={Wall-clock time [s]},\n');
    fprintf (fid, '    title={Cost vs max\\_level},\n');
    fprintf (fid, '    legend pos=outer north east, grid=both, xtick=%s]\n', xtick_str);
    for im = 1:numel(methods)
        fprintf (fid, '\\addplot[%s, mark=%s, mark size=2.5pt, very thick] table[x=max_level, y=wallclock] {adaptivity_summary_%s.dat};\n', ...
                 method_colors{1+mod(im-1,numel(method_colors))}, ...
                 method_marks {1+mod(im-1,numel(method_marks))}, methods{im});
        fprintf (fid, '\\addlegendentry{%s}\n', method_legend(methods{im}));
    end
    fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');

    % ---- Load-displacement: one tikzpicture per method, all levels overlaid
    for im = 1:numel(methods)
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{axis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={$u_r$ at outer boundary [mm]},\n');
        fprintf (fid, '    ylabel={Internal pressure $P$ [MPa]},\n');
        fprintf (fid, '    title={Load--displacement, %s},\n', method_legend(methods{im}));
        fprintf (fid, '    legend pos=outer north east, grid=both]\n');
        fprintf (fid, '\\addplot[black, thick, no marks] table[x=u, y=P] {adaptivity_load_disp_hill.dat};\n');
        fprintf (fid, '\\addlegendentry{Hill}\n');
        for il = 1:numel(max_levels)
            ml = max_levels(il);
            col = level_colors{1+mod(il-1,numel(level_colors))};
            mk  = level_marks {1+mod(il-1,numel(level_marks))};
            fprintf (fid, '\\addplot[%s, mark=%s, mark size=2pt, very thick] table[x=u, y=P] {adaptivity_load_disp_%s_lvl%d.dat};\n', ...
                     col, mk, methods{im}, ml);
            fprintf (fid, '\\addlegendentry{max\\_level = %d}\n', ml);
        end
        fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    % ---- Stress profiles: one row of two tikzpictures per method --------
    for im = 1:numel(methods)
        for kind = 1:2  % 1 = sigma_r, 2 = sigma_t
            if kind == 1
                ylab = '$\sigma_r$ [MPa]'; ttl = 'Radial';
                hill_file = sprintf('adaptivity_sigma_r_hill_step%d.dat', last_step);
                fname_pat = 'adaptivity_sigma_r_%s_lvl%d.dat';
            else
                ylab = '$\sigma_t$ [MPa]'; ttl = 'Tangential';
                hill_file = sprintf('adaptivity_sigma_t_hill_step%d.dat', last_step);
                fname_pat = 'adaptivity_sigma_t_%s_lvl%d.dat';
            end
            fprintf (fid, '\\begin{tikzpicture}\n');
            fprintf (fid, '\\begin{axis}[\n');
            fprintf (fid, '    width=10cm, height=7cm,\n');
            fprintf (fid, '    xlabel={Radial coordinate $r$ [mm]},\n');
            fprintf (fid, '    ylabel={%s},\n', ylab);
            fprintf (fid, '    title={%s stress, %s, step %d},\n', ttl, method_legend(methods{im}), last_step);
            fprintf (fid, '    legend pos=outer north east, grid=both]\n');
            fprintf (fid, '\\addplot[black, thick, no marks] table[x=r, y=hill] {%s};\n', hill_file);
            fprintf (fid, '\\addlegendentry{Hill}\n');
            for il = 1:numel(max_levels)
                ml = max_levels(il);
                col = level_colors{1+mod(il-1,numel(level_colors))};
                mk  = level_marks {1+mod(il-1,numel(level_marks))};
                fprintf (fid, '\\addplot[%s, mark=%s, mark size=2pt, mark repeat=10, very thick] table[x=r, y=sigma] {%s};\n', ...
                         col, mk, sprintf(fname_pat, methods{im}, ml));
                fprintf (fid, '\\addlegendentry{max\\_level = %d}\n', ml);
            end
            fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
        end
    end

    % ---- Comparison: all methods at highest level on one plot -----------
    best_lvl = max(max_levels);
    for kind = 1:2
        if kind == 1
            ylab = '$\sigma_r$ [MPa]'; ttl = 'Radial';
            hill_file = sprintf('adaptivity_sigma_r_hill_step%d.dat', last_step);
            fname_pat = 'adaptivity_sigma_r_%s_lvl%d.dat';
        else
            ylab = '$\sigma_t$ [MPa]'; ttl = 'Tangential';
            hill_file = sprintf('adaptivity_sigma_t_hill_step%d.dat', last_step);
            fname_pat = 'adaptivity_sigma_t_%s_lvl%d.dat';
        end
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{axis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={Radial coordinate $r$ [mm]},\n');
        fprintf (fid, '    ylabel={%s},\n', ylab);
        fprintf (fid, '    title={%s stress, step %d},\n', ttl, last_step);
        fprintf (fid, '    legend pos=outer north east, grid=both]\n');
        fprintf (fid, '\\addplot[black, thick, no marks] table[x=r, y=hill] {%s};\n', hill_file);
        fprintf (fid, '\\addlegendentry{Hill}\n');
        for im = 1:numel(methods)
            fprintf (fid, '\\addplot[%s, mark=%s, mark size=2pt, mark repeat=10, very thick] table[x=r, y=sigma] {%s};\n', ...
                     method_colors{1+mod(im-1,numel(method_colors))}, ...
                     method_marks {1+mod(im-1,numel(method_marks))}, ...
                     sprintf(fname_pat, methods{im}, best_lvl));
            fprintf (fid, '\\addlegendentry{%s}\n', method_legend(methods{im}));
        end
        fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    fprintf (fid, '\\end{document}\n');
end

function s = latex_escape(s)
    s = strrep(s, '_', '\_');
end

function plot_adaptivity_results (R, methods, max_levels, u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius, last_step)
    figure(1); clf;
    cmap = lines(numel(methods));
    metrics = {'err_u_inf','err_u_avg','err_sigma_r_inf','err_sigma_r_avg','err_sigma_t_inf','err_sigma_t_avg'};
    titles  = {'max u-error','avg u-error','max \sigma_r error','avg \sigma_r error','max \sigma_t error','avg \sigma_t error'};
    for k = 1:numel(metrics)
        subplot(2,3,k); hold on; grid on;
        for im = 1:numel(methods)
            semilogy(max_levels, R.(metrics{k})(im,:), '-o', 'Color', cmap(im,:), ...
                     'LineWidth', 1.4, 'DisplayName', methods{im});
        end
        set(gca, 'YScale', 'log');
        xlabel('max\_level'); ylabel(titles{k});
        legend('Location', 'best');
    end
    saveas(figure(1), 'adaptivity_errors.png');
end
