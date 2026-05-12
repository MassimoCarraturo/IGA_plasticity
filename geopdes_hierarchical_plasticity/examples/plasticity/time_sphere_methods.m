% TIME_SPHERE_METHODS  Measure wall-clock and CPU time for each projection
% method (L2, QI, QI_C0), both with the libqi backend enabled and with it
% bypassed (MATLAB reference for the local-LS step).
%
% Writes sphere_timings.mat and sphere_timings.csv.

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');

% Path setup
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

% --- problem setup (same as ex_plastic_sphere_hier) ----------------------
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

p = 2;
method_data_template.degree     = [p p p];
method_data_template.regularity = [p-1 p-1 p-1];
method_data_template.nsub_coarse = [15, 5, 5];
method_data_template.nsub_refine = [2 2 2];
method_data_template.nquad      = [p+1 p+1 p+1];
method_data_template.space_type = 'standard';
method_data_template.truncated  = 1;
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

% Locate the libqi MEX file so we can hide / restore it to switch backends
mex_dir  = fullfile(project_root, 'geopdes_hierarchical_plasticity', ...
                    'quasi_interpolant_hierarchical', 'libqi', 'matlab');
mex_file = fullfile(mex_dir, ['qi_local_ls_mex.', mexext()]);
mex_off  = [mex_file, '.disabled'];
if ~exist(mex_file, 'file')
    if exist(mex_off, 'file')
        movefile(mex_off, mex_file);
    else
        error('time_sphere_methods:mex', ...
              'qi_local_ls_mex MEX file missing — run libqi/build_mex.m first');
    end
end

% Verify both modes work by toggling once:
fprintf('libqi MEX path: %s\n', mex_file);

methods   = {'L2', 'QI', 'QI_C0'};
backends  = {'matlab', 'libqi'};
results   = struct();

for ib = 1:numel(backends)
    backend = backends{ib};
    if strcmp(backend, 'matlab') && exist(mex_file, 'file')
        movefile(mex_file, mex_off);   % hide the MEX
    elseif strcmp(backend, 'libqi') && ~exist(mex_file, 'file')
        movefile(mex_off, mex_file);   % restore
    end

    % Force the persistent dispatch in history_variable_projection_hier to be
    % re-resolved on next call. A fresh `clear functions` invalidates persistents.
    clear functions; %#ok<CLFUNC>
    addpath (fullfile(project_root, 'geopdes_hierarchical_plasticity', ...
             'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

    if exist('qi_local_ls_mex', 'file') == 3
        fprintf('\n==== backend = libqi (MEX present) ====\n');
    else
        fprintf('\n==== backend = matlab (MEX hidden) ====\n');
    end

    for im = 1:numel(methods)
        method = methods{im};
        fprintf('--- %-7s with %s ---\n', method, backend);

        method_data = method_data_template;
        method_data.type_projection = method;

        cpu0  = cputime();
        wall0 = tic;
        [~, ~, ~, ~, ~, ~, ~, ~] = adaptivity_J2_plasticity ...
            (problem_data, method_data, adaptivity_data);
        wall  = toc(wall0);
        cpu   = cputime() - cpu0;

        results.(matlab.lang.makeValidName([method '_' backend])) = ...
            struct('method', method, 'backend', backend, ...
                   'wallclock_s', wall, 'cpu_s', cpu);
        fprintf('    wallclock = %.2f s   CPU = %.2f s   parallel-eff = %.2fx\n', ...
                wall, cpu, cpu / max(wall, eps));
    end
end

% Restore MEX if it was hidden
if exist(mex_off, 'file')
    movefile(mex_off, mex_file);
end

save (fullfile(here, 'sphere_timings.mat'), '-struct', 'results');

% --- write CSV table -----------------------------------------------------
fid = fopen (fullfile(here, 'sphere_timings.csv'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'method,backend,wallclock_s,cpu_s,cpu_per_wall\n');
fields = fieldnames(results);
for i = 1:numel(fields)
    R = results.(fields{i});
    fprintf(fid, '%s,%s,%.3f,%.3f,%.2f\n', R.method, R.backend, ...
            R.wallclock_s, R.cpu_s, R.cpu_s / max(R.wallclock_s, eps));
end
fprintf('Wrote sphere_timings.mat and sphere_timings.csv\n');
