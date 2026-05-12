% COMPARE_ESTIMATORS  Compare residual-divergence and stress-gradient
% error estimators at fixed max_level = 3 across the three projection
% methods (L2 / QI / QI_C0).  Uses the same coarser-mesh setup
% (nsub_coarse = [3 1 1], num_max_iter = 5) as adaptivity_max_level_study.
%
% Reports max & average errors against Hill, plus runtime, for each
% (estimator, method) combination.  Writes
%   compare_estimators.csv
%   compare_estimators.mat
%   pgfplots/compare_estimators.dat
%   pgfplots/compare_estimators.tex

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
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

% Problem data (same as the rest of the study)
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
nload = 5; P = 332*.99;
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.p = @(x, y, z) P*ones (size (x));

[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution (E, nu, problem_data.yield_stress(1), 100, 200, P, nload);
load_step_eval_stress = [2, 5];

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

adaptivity_data_template.flag = 'elements';
adaptivity_data_template.C0_est = 1.0;
adaptivity_data_template.mark_param = 0.9;
adaptivity_data_template.mark_param_coarsening = 0.05;
adaptivity_data_template.mark_strategy = 'MS';
adaptivity_data_template.max_ndof = 200000;
adaptivity_data_template.max_nel  = 100000;
adaptivity_data_template.num_max_iter = 5;
adaptivity_data_template.tol = 1e-5 * 3.14 * 30000;
adaptivity_data_template.adm_strategy = 'admissible';
adaptivity_data_template.coarsening_flag = 'any';
adaptivity_data_template.adm = p;
adaptivity_data_template.max_level = 3;   % fixed for the comparison

methods    = {'L2', 'QI', 'QI_C0'};
estimators = {'div_sigma', 'stress_gradient'};

R = struct();
R.methods    = methods;
R.estimators = estimators;
fields = {'wall','cpu','ndof','nel','nlevels','err_u_max','err_u_avg', ...
          'err_sigma_r_max','err_sigma_r_avg','err_sigma_t_max','err_sigma_t_avg'};
for k = 1:numel(fields)
    R.(fields{k}) = nan(numel(methods), numel(estimators));
end

for ie = 1:numel(estimators)
    for im = 1:numel(methods)
        fprintf('\n========================================================\n');
        fprintf('  estimator = %-16s   method = %s   (%d/%d)\n', ...
                estimators{ie}, methods{im}, ...
                (ie-1)*numel(methods)+im, numel(methods)*numel(estimators));
        fprintf('========================================================\n');

        method_data = method_data_template;
        method_data.type_projection = methods{im};
        adaptivity_data = adaptivity_data_template;
        adaptivity_data.estimator = estimators{ie};

        cpu0 = cputime(); wall0 = tic;
        try
            [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, cell_u, ~, cell_sigma, ~] = ...
                adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);
        catch e
            fprintf('  FAILED: %s\n', e.message);
            R.wall(im, ie) = toc(wall0);
            R.cpu(im, ie)  = cputime() - cpu0;
            save (fullfile(here, 'compare_estimators.mat'), '-struct', 'R');
            continue;
        end
        R.wall(im, ie) = toc(wall0);
        R.cpu(im, ie)  = cputime() - cpu0;

        R.ndof(im, ie)    = cell_hspace{end}.ndof;
        R.nel(im, ie)     = cell_hmsh{end}.nel;
        R.nlevels(im, ie) = cell_hmsh{end}.nlevels;

        u_r = zeros(1, nload+1);
        P_i = zeros(1, nload+1);
        for i = 1:nload+1
            P_i(i) = P/nload * (i-1);
            if i > 1
                eu = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1, .5, .5});
                u_r(i) = norm(eu);
            end
        end
        keep = [true, diff(P_ex) > 0];
        Pi_pos = P_i(P_i > 0);  u_r_pos = u_r(P_i > 0);
        hill_u = interp1(P_ex(keep), u_ex(keep), Pi_pos, 'linear', 'extrap');
        diff_u = abs(u_r_pos - hill_u);
        R.err_u_max(im, ie) = max(diff_u);
        R.err_u_avg(im, ie) = mean(diff_u);

        rad_dir = [.5; .5; sqrt(2)/2];
        tan_dir = [.5; .5; -sqrt(2)/2];
        voigt   = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];
        rad_pos = linspace(0, 1, length(radius));

        sigma_rad = zeros(length(radius), length(load_step_eval_stress));
        sigma_tan = zeros(length(radius), length(load_step_eval_stress));
        for jload = 1:length(load_step_eval_stress)
            for jpoint = 1:length(rad_pos)
                position = {rad_pos(jpoint), 0.5, 0.75};
                M = zeros(3,3);
                for icomp = 1:6
                    eu = sp_eval (cell_sigma{load_step_eval_stress(jload)}(:,icomp), ...
                                  cell_hspace_scalar{load_step_eval_stress(jload)}, geometry, position);
                    M(voigt(icomp,1), voigt(icomp,2)) = eu;
                    if icomp > 3, M(voigt(icomp,2), voigt(icomp,1)) = eu; end
                end
                sigma_rad(jpoint, jload) = rad_dir' * M * rad_dir;
                sigma_tan(jpoint, jload) = tan_dir' * M * tan_dir;
            end
        end
        diff_sr = zeros(length(radius), length(load_step_eval_stress));
        diff_st = zeros(length(radius), length(load_step_eval_stress));
        for jload = 1:length(load_step_eval_stress)
            ls = load_step_eval_stress(jload);
            diff_sr(:,jload) = abs(sigma_rad(:,jload) - sigma_r_ex(ls,:).');
            diff_st(:,jload) = abs(sigma_tan(:,jload) - sigma_t_ex(ls,:).');
        end
        R.err_sigma_r_max(im, ie) = max(diff_sr(:));
        R.err_sigma_r_avg(im, ie) = mean(diff_sr(:));
        R.err_sigma_t_max(im, ie) = max(diff_st(:));
        R.err_sigma_t_avg(im, ie) = mean(diff_st(:));

        fprintf('  ndof=%d  nel=%d  nlev=%d  wall=%.1fs\n', ...
                R.ndof(im, ie), R.nel(im, ie), R.nlevels(im, ie), R.wall(im, ie));
        fprintf('  err_u   max=%.3e  avg=%.3e\n', R.err_u_max(im,ie), R.err_u_avg(im,ie));
        fprintf('  sigma_r max=%.3e  avg=%.3e\n', R.err_sigma_r_max(im,ie), R.err_sigma_r_avg(im,ie));
        fprintf('  sigma_t max=%.3e  avg=%.3e\n', R.err_sigma_t_max(im,ie), R.err_sigma_t_avg(im,ie));

        save (fullfile(here, 'compare_estimators.mat'), '-struct', 'R');
    end
end

% CSV summary
fid = fopen (fullfile(here, 'compare_estimators.csv'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf (fid, 'estimator,method,ndof,nel,nlevels,wall_s,cpu_s,err_u_max,err_u_avg,err_sr_max,err_sr_avg,err_st_max,err_st_avg\n');
for ie = 1:numel(estimators)
    for im = 1:numel(methods)
        fprintf (fid, '%s,%s,%d,%d,%d,%.3f,%.3f,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n', ...
                 estimators{ie}, methods{im}, ...
                 R.ndof(im, ie), R.nel(im, ie), R.nlevels(im, ie), ...
                 R.wall(im, ie), R.cpu(im, ie), ...
                 R.err_u_max(im,ie), R.err_u_avg(im,ie), ...
                 R.err_sigma_r_max(im,ie), R.err_sigma_r_avg(im,ie), ...
                 R.err_sigma_t_max(im,ie), R.err_sigma_t_avg(im,ie));
    end
end
clear cleanup;

% pgfplots: a bar-chart-like table per estimator
out_dir = 'pgfplots';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
for ie = 1:numel(estimators)
    method_idx = (1:numel(methods)).';
    rows = [method_idx, R.err_u_max(:,ie), R.err_u_avg(:,ie), ...
            R.err_sigma_r_max(:,ie), R.err_sigma_r_avg(:,ie), ...
            R.err_sigma_t_max(:,ie), R.err_sigma_t_avg(:,ie), ...
            R.wall(:,ie), R.ndof(:,ie), R.nel(:,ie)];
    write_dat (fullfile(out_dir, sprintf('compare_%s.dat', estimators{ie})), ...
               {'method_idx','err_u_max','err_u_avg', ...
                'err_sr_max','err_sr_avg','err_st_max','err_st_avg', ...
                'wall','ndof','nel'}, rows);
end

write_pgfplots (fullfile(out_dir, 'compare_estimators.tex'), methods, estimators);

fprintf('\nWrote compare_estimators.csv / .mat and %s\n', ...
        fullfile(out_dir, 'compare_estimators.tex'));

% =========================================================================
function write_dat (filename, headers, data)
    fid = fopen(filename, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', headers{1});
    for k = 2:numel(headers); fprintf(fid, '\t%s', headers{k}); end
    fprintf(fid, '\n');
    % '\t' / '\n' in single-quoted MATLAB are 2-char literals; build the
    % format string explicitly so fprintf interprets them as tab/newline.
    fmt = [repmat('%.10g\t', 1, size(data, 2) - 1), '%.10g\n'];
    for i = 1:size(data,1); fprintf(fid, fmt, data(i,:)); end
end

function write_pgfplots (filename, methods, estimators)
    fid = fopen(filename, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%% Auto-generated by compare_estimators.m\n');
    fprintf(fid, '\\documentclass[tikz,border=4pt]{standalone}\n');
    fprintf(fid, '\\usepackage{pgfplots}\n');
    fprintf(fid, '\\pgfplotsset{compat=1.18, every axis/.append style={font=\\small}}\n\n');
    fprintf(fid, '\\begin{document}\n\n');

    % NOTE: single-quoted MATLAB strings do NOT process backslash escapes,
    % so '\mathrm' is the 7-char literal that fprintf writes verbatim.
    % Keep them with single backslashes — '\\' here would write '\\' to
    % the .tex (which LaTeX parses as a tabular line-break, not as a
    % command).
    metrics = {
        'err_u_max',  'Max $|u_r - u_\mathrm{Hill}|$',     'Max u-error';
        'err_u_avg',  'Avg $|u_r - u_\mathrm{Hill}|$',     'Avg u-error';
        'err_sr_max', 'Max $|\sigma_r - \sigma_r^\mathrm{Hill}|$ [MPa]', 'Max $\sigma_r$ error';
        'err_sr_avg', 'Avg $|\sigma_r - \sigma_r^\mathrm{Hill}|$ [MPa]', 'Avg $\sigma_r$ error';
        'err_st_max', 'Max $|\sigma_t - \sigma_t^\mathrm{Hill}|$ [MPa]', 'Max $\sigma_t$ error';
        'err_st_avg', 'Avg $|\sigma_t - \sigma_t^\mathrm{Hill}|$ [MPa]', 'Avg $\sigma_t$ error';
    };

    est_colors = {'blue!60!black', 'orange!80!black'};
    est_marks  = {'*', 'square*'};
    method_labels = cell(numel(methods), 1);
    for i = 1:numel(methods)
        method_labels{i} = sprintf('{%s}', strrep(methods{i}, '_', '\_'));
    end
    xticklabels_list = strjoin(method_labels, ',');

    for k = 1:size(metrics,1)
        col   = metrics{k,1};
        ylab  = metrics{k,2};
        ttl   = metrics{k,3};

        fprintf(fid, '\\begin{tikzpicture}\n');
        fprintf(fid, '\\begin{axis}[\n');
        fprintf(fid, '  ybar, ymode=log, log origin=infty,\n');
        fprintf(fid, '  width=11cm, height=7cm, bar width=12pt,\n');
        fprintf(fid, '  xtick={1,2,3}, xticklabels={%s},\n', xticklabels_list);
        fprintf(fid, '  enlarge x limits=0.2,\n');
        fprintf(fid, '  ylabel={%s},\n', ylab);
        fprintf(fid, '  title={%s vs Hill, max\\_level = 3},\n', ttl);
        fprintf(fid, '  legend pos=outer north east, grid=both]\n');
        for ie = 1:numel(estimators)
            fprintf(fid, '\\addplot+[fill=%s, draw=%s] table[x=method_idx, y=%s] {compare_%s.dat};\n', ...
                    est_colors{1+mod(ie-1,numel(est_colors))}, ...
                    est_colors{1+mod(ie-1,numel(est_colors))}, ...
                    col, estimators{ie});
            fprintf(fid, '\\addlegendentry{%s}\n', strrep(estimators{ie}, '_', '\_'));
        end
        fprintf(fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    % Wall-clock per estimator (linear axis)
    fprintf(fid, '\\begin{tikzpicture}\n');
    fprintf(fid, '\\begin{axis}[\n');
    fprintf(fid, '  ybar, width=11cm, height=7cm, bar width=12pt,\n');
    fprintf(fid, '  xtick={1,2,3}, xticklabels={%s},\n', xticklabels_list);
    fprintf(fid, '  enlarge x limits=0.2,\n');
    fprintf(fid, '  ylabel={Wall-clock [s]},\n');
    fprintf(fid, '  title={Cost per (estimator, method), max\\_level = 3},\n');
    fprintf(fid, '  legend pos=outer north east, grid=both]\n');
    for ie = 1:numel(estimators)
        fprintf(fid, '\\addplot+[fill=%s, draw=%s] table[x=method_idx, y=wall] {compare_%s.dat};\n', ...
                est_colors{1+mod(ie-1,numel(est_colors))}, ...
                est_colors{1+mod(ie-1,numel(est_colors))}, ...
                estimators{ie});
        fprintf(fid, '\\addlegendentry{%s}\n', strrep(estimators{ie}, '_', '\_'));
    end
    fprintf(fid, '\\end{axis}\n\\end{tikzpicture}\n\n');

    fprintf(fid, '\\end{document}\n');
end
