% BENCHMARK_OLD_VS_NEW  Head-to-head timing: old libqi v1.0 vs optimised v1.1.
%
%   Swaps between qi_mex_old.mexw64 / qi_local_ls_mex_old.mexw64  (baseline)
%   and the optimised .mexw64 builds, running the same sphere benchmark with
%   each.  Only QI and QI_C0 methods are timed (L2 doesn't use libqi).
%
%   Output: console table + benchmark_old_vs_new.mat

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

mex_dir = fullfile(project_root, 'geopdes_hierarchical_plasticity', ...
    'quasi_interpolant_hierarchical', 'libqi', 'matlab');

% ========================================================================
% Problem setup (same as benchmark_qi_optimized.m)
% ========================================================================
problem_data.geo_name = 'geo_eighth_sphere.txt';
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [1];
problem_data.symm_sides   = [3 4 5];
problem_data.slider_sides = [6];
problem_data.penalty_slider = @(x, y, z) 1e8 * ones(size(x));

E  = 210000; nu = 0.3;
problem_data.yield_stress = @(x, y, z) 240 * ones(size(x));
problem_data.kappa_lame   = @(x, y, z) E/(3*(1-2*nu)) * ones(size(x));
problem_data.mu_lame      = @(x, y, z) E/(2*(1+nu)) * ones(size(x));

nload = 3;
P     = 332 * 0.99;
problem_data.f = @(x, y, z) zeros(3, size(x,1), size(x,2), size(x,3));
problem_data.g = @(x, y, z, ind) zeros(3, size(x,1), size(x,2), size(x,3));
problem_data.h = @(x, y, z, ind) zeros(3, size(x,1), size(x,2), size(x,3));
problem_data.p = @(x, y, z) P * ones(size(x));

p = 2;
method_data_template.degree      = [p p p];
method_data_template.regularity  = [p-1 p-1 p-1];
method_data_template.nsub_coarse = [4, 4, 4];
method_data_template.nsub_refine = [2 2 2];
method_data_template.nquad       = [p+1 p+1 p+1];
method_data_template.space_type  = 'standard';
method_data_template.truncated   = 1;
method_data_template.nload       = nload;
method_data_template.newton_tol  = 1e-6;
method_data_template.newton_tol_abs = 1e-8;
method_data_template.newton_iter_max = 50;

adaptivity_data.flag        = 'elements';
adaptivity_data.estimator   = 'sphere_front';
adaptivity_data.C0_est      = 1.0;
adaptivity_data.mark_param  = 0.9;
adaptivity_data.mark_param_coarsening = 0.05;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level   = 3;
adaptivity_data.max_ndof    = 10000;
adaptivity_data.num_max_iter = 2;
adaptivity_data.max_nel     = 4000;
adaptivity_data.tol         = 1e-5 * 3.14 * 30000;
adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any';
adaptivity_data.adm         = p;

% ========================================================================
% Helper: swap MEX files
% ========================================================================
function swap_mex(mex_dir, tag)
% tag = 'old' or 'new'
    qi_new  = fullfile(mex_dir, 'qi_mex.mexw64');
    qi_old  = fullfile(mex_dir, 'qi_mex_old.mexw64');
    qi_save = fullfile(mex_dir, 'qi_mex_new.mexw64');
    ls_new  = fullfile(mex_dir, 'qi_local_ls_mex.mexw64');
    ls_old  = fullfile(mex_dir, 'qi_local_ls_mex_old.mexw64');
    ls_save = fullfile(mex_dir, 'qi_local_ls_mex_new.mexw64');

    % Clear MEX from memory so MATLAB reloads the swapped binary
    clear qi_mex qi_local_ls_mex;

    if strcmpi(tag, 'old')
        % Save current (new) as *_new, copy old → active
        if isfile(qi_save);  delete(qi_save);  end
        if isfile(ls_save);  delete(ls_save);  end
        copyfile(qi_new, qi_save);
        copyfile(ls_new, ls_save);
        copyfile(qi_old, qi_new);
        copyfile(ls_old, ls_new);
    else  % 'new'
        % Restore new → active
        if isfile(qi_save)
            copyfile(qi_save, qi_new);
            copyfile(ls_save, ls_new);
            delete(qi_save);
            delete(ls_save);
        end
    end
end

% ========================================================================
% Run benchmark
% ========================================================================
proj_methods = {'QI', 'QI_C0'};
versions     = {'old', 'new'};
results = struct();

for iv = 1:numel(versions)
    ver = versions{iv};
    swap_mex(mex_dir, ver);

    for imethod = 1:numel(proj_methods)
        method_name = proj_methods{imethod};
        label = sprintf('%s_%s', method_name, ver);
        fprintf('\n================================================================\n');
        fprintf('  %s : type_projection = %s, libqi = %s\n', label, method_name, ver);
        fprintf('================================================================\n');

        method_data = method_data_template;
        method_data.type_projection = method_name;

        % Clear persistent handles so the newly-swapped MEX is picked up
        clear history_variable_projection_hier;

        t0 = tic;
        [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
         cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
            adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);
        t_wall = toc(t0);

        ndof_final = cell_hspace{nload}.ndof;
        nel_final  = sum(cell_hmsh{nload}.nel_per_level);

        R = struct();
        R.method    = method_name;
        R.version   = ver;
        R.wall_time = t_wall;
        R.ndof      = ndof_final;
        R.nel       = nel_final;
        results.(matlab.lang.makeValidName(label)) = R;

        fprintf('  %s done: wall = %.2f s, ndof = %d, nel = %d\n', ...
                label, t_wall, ndof_final, nel_final);
    end
end

% Ensure we end with the NEW (optimised) version active
swap_mex(mex_dir, 'new');

% ========================================================================
% Summary
% ========================================================================
fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  OLD vs NEW libqi — head-to-head comparison\n');
fprintf('================================================================\n');
fprintf('%-12s  %10s  %8s  %8s\n', 'Run', 'WallTime', 'ndof', 'nel');
fprintf('%-12s  %10s  %8s  %8s\n', '---', '--------', '----', '---');

fields = fieldnames(results);
for i = 1:numel(fields)
    R = results.(fields{i});
    fprintf('%-12s  %9.2f s  %8d  %8d\n', fields{i}, R.wall_time, R.ndof, R.nel);
end

fprintf('\nSpeedup (new / old):\n');
for imethod = 1:numel(proj_methods)
    mn = proj_methods{imethod};
    old_label = matlab.lang.makeValidName(sprintf('%s_old', mn));
    new_label = matlab.lang.makeValidName(sprintf('%s_new', mn));
    if isfield(results, old_label) && isfield(results, new_label)
        t_old = results.(old_label).wall_time;
        t_new = results.(new_label).wall_time;
        fprintf('  %-8s: old = %.2f s, new = %.2f s  =>  %.2fx speedup\n', ...
                mn, t_old, t_new, t_old / t_new);
    end
end
fprintf('\n');

save(fullfile(results_dir, 'benchmark_old_vs_new.mat'), 'results', 'proj_methods', 'versions');
fprintf('Results saved to %s\n', fullfile(results_dir, 'benchmark_old_vs_new.mat'));

% ========================================================================
% Local function: swap MEX files
% ========================================================================
