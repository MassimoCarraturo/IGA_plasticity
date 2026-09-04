function run_ci_tests(suite)
% RUN_CI_TESTS  CI test runner — builds MEX, then runs selected test suites.
%
%   run_ci_tests()           — runs all suites: build, unit, smoke
%   run_ci_tests('build')    — build MEX only
%   run_ci_tests('unit')     — unit tests only  (test_localLS_vs_matlab)
%   run_ci_tests('smoke')    — smoke integration (smoke_test_estimator)
%   run_ci_tests('benchmark')— performance tests (benchmark_qi_optimized)
%   run_ci_tests('all')      — everything including benchmark
%
%   Exit codes:  0 = all passed,  1 = one or more failures.
%   (Called via `matlab -batch "run_ci_tests"` so MATLAB exits with the code.)

    if nargin < 1, suite = 'default'; end

    project_root = fileparts(fileparts(mfilename('fullpath')));
    setup_paths(project_root);

    results = {};
    t_total = tic;

    switch lower(suite)
        case 'build'
            results{end+1} = run_build(project_root);
        case 'unit'
            results{end+1} = run_unit_tests(project_root);
        case 'smoke'
            results{end+1} = run_smoke_tests(project_root);
        case 'benchmark'
            results{end+1} = run_benchmark(project_root);
        case 'all'
            results{end+1} = run_build(project_root);
            results{end+1} = run_unit_tests(project_root);
            results{end+1} = run_smoke_tests(project_root);
            results{end+1} = run_benchmark(project_root);
        otherwise  % 'default' — build + unit + smoke
            results{end+1} = run_build(project_root);
            results{end+1} = run_unit_tests(project_root);
            results{end+1} = run_smoke_tests(project_root);
    end

    % ---- Summary ----
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  CI TEST SUMMARY\n');
    fprintf('================================================================\n');
    n_pass = 0; n_fail = 0;
    for k = 1:numel(results)
        R = results{k};
        if R.passed
            status = 'PASS';
            n_pass = n_pass + 1;
        else
            status = 'FAIL';
            n_fail = n_fail + 1;
        end
        fprintf('  [%s]  %-40s  (%.1f s)\n', status, R.name, R.time);
        if ~isempty(R.message)
            fprintf('         %s\n', R.message);
        end
    end
    fprintf('----------------------------------------------------------------\n');
    fprintf('  Total: %d passed, %d failed  (%.1f s)\n', ...
            n_pass, n_fail, toc(t_total));
    fprintf('================================================================\n');

    if n_fail > 0
        error('run_ci_tests:failed', '%d test suite(s) failed.', n_fail);
    end
end

% ========================================================================
%  Path setup
% ========================================================================
function setup_paths(project_root)
    addpath(genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
    addpath(genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
    addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
    addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
    addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'initialize')));
    addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
    addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
    addpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', ...
        'quasi_interpolant_hierarchical', 'libqi', 'matlab'));
end

% ========================================================================
%  Suite: Build MEX
% ========================================================================
function R = run_build(project_root)
    R.name = 'Build MEX (qi_mex + qi_local_ls_mex)';
    t0 = tic;
    try
        old_dir = pwd;
        cleanup = onCleanup(@() cd(old_dir));
        build_mex_path = fullfile(project_root, ...
            'geopdes_hierarchical_plasticity', ...
            'quasi_interpolant_hierarchical', 'libqi');
        cd(build_mex_path);
        build_mex();
        % Verify the MEX files exist
        mex_dir = fullfile(build_mex_path, 'matlab');
        ext = mexext();
        assert(isfile(fullfile(mex_dir, ['qi_mex.', ext])), ...
               'qi_mex.%s not found after build', ext);
        assert(isfile(fullfile(mex_dir, ['qi_local_ls_mex.', ext])), ...
               'qi_local_ls_mex.%s not found after build', ext);
        R.passed  = true;
        R.message = sprintf('Built qi_mex.%s and qi_local_ls_mex.%s', ext, ext);
    catch ME
        R.passed  = false;
        R.message = ME.message;
    end
    R.time = toc(t0);
end

% ========================================================================
%  Suite: Unit tests (test_localLS_vs_matlab)
% ========================================================================
function R = run_unit_tests(project_root)
    R.name = 'Unit: libqi MEX vs MATLAB reference';
    t0 = tic;
    try
        warning('off', 'MATLAB:nearlySingularMatrix');
        warning('off', 'MATLAB:singularMatrix');

        [n_pass, n_fail, details] = test_localLS_ci(project_root);

        if n_fail == 0
            R.passed  = true;
            R.message = sprintf('%d/%d checks passed', n_pass, n_pass + n_fail);
        else
            R.passed  = false;
            R.message = sprintf('%d/%d checks FAILED: %s', ...
                        n_fail, n_pass + n_fail, strjoin(details, '; '));
        end
    catch ME
        R.passed  = false;
        R.message = ME.message;
    end
    R.time = toc(t0);
end

function [n_pass, n_fail, details] = test_localLS_ci(project_root)
% Reimplements test_localLS_vs_matlab logic but returns structured results.
    n_pass = 0; n_fail = 0; details = {};

    % Build spaces
    deg  = [3 3];
    nelx = 8; nely = 8;
    geometry = geo_load(nrb4surf([0 0], [1 0], [0 1], [1 1]));
    [knots, zeta] = kntrefine(geometry.nurbs.knots, [nelx nely]-1, deg, deg-1);
    rule = msh_gauss_nodes(deg + 1);
    [qn, qw] = msh_set_quad_nodes(zeta, rule);
    msh_obj  = msh_cartesian(zeta, qn, qw, geometry);
    space    = sp_bspline(knots, deg, msh_obj);
    hmsh     = hierarchical_mesh(msh_obj, [2 2]);
    hspace   = hierarchical_space(hmsh, space, 'standard', false, deg-1);

    % Regularity-0 space
    [knots0, zeta0] = kntrefine(geometry.nurbs.knots, [nelx nely]-1, deg, [0 0]);
    [qn0, qw0] = msh_set_quad_nodes(zeta0, rule);
    msh0   = msh_cartesian(zeta0, qn0, qw0, geometry);
    space0 = sp_bspline(knots0, deg, msh0);
    hmsh0  = hierarchical_mesh(msh0, [2 2]);
    hspace0 = hierarchical_space(hmsh0, space0, 'standard', false, [0 0]);

    % Generate data
    rng(0);
    [gx, gy] = ndgrid(linspace(0.01, 0.99, 50), linspace(0.01, 0.99, 50));
    x = gx(:) + 0.005*(rand(numel(gx),1)-0.5);
    y = gy(:) + 0.005*(rand(numel(gy),1)-0.5);
    fv = [(tanh(9*y-9*x)+1)/9, exp(-((10*x-6).^2+(10*y+7).^2))/1.5, ...
          sin(3*pi*x).*cos(2*pi*y)];
    data = [x, y];

    tol = 1e-8;

    % Standard regularity tests
    for lambda = [0, 1e-9, 1e-3]
        label = sprintf('reg=default, lambda=%g', lambda);
        QI_ref = getcoeff_localLS_Bspl(hspace, hmsh, data, fv, lambda);
        QI_c   = getcoeff_localLS_Bspl_c(hspace, hmsh, data, fv, lambda);
        rel = max(abs(QI_ref(:)-QI_c(:))) / max(abs(QI_ref(:)));
        fprintf('  [%s] rel_err = %.3e\n', label, rel);
        if rel < tol
            n_pass = n_pass + 1;
        else
            n_fail = n_fail + 1;
            details{end+1} = sprintf('%s: rel=%.2e', label, rel); %#ok<AGROW>
        end
    end

    % Regularity-0 tests
    for lambda = [0, 1e-9, 1e-3]
        label = sprintf('reg=0, lambda=%g', lambda);
        QI_ref = getcoeff_localLS_Bspl(hspace0, hmsh0, data, fv, lambda);
        QI_c   = getcoeff_localLS_Bspl_c(hspace0, hmsh0, data, fv, lambda);
        rel = max(abs(QI_ref(:)-QI_c(:))) / max(abs(QI_ref(:)));
        fprintf('  [%s] rel_err = %.3e\n', label, rel);
        if rel < tol
            n_pass = n_pass + 1;
        else
            n_fail = n_fail + 1;
            details{end+1} = sprintf('%s: rel=%.2e', label, rel); %#ok<AGROW>
        end
    end
end

% ========================================================================
%  Suite: Smoke integration test
% ========================================================================
function R = run_smoke_tests(project_root)
    R.name = 'Smoke: adaptive J2 plasticity (sphere, QI)';
    t0 = tic;
    try
        warning('off', 'MATLAB:nearlySingularMatrix');
        warning('off', 'MATLAB:singularMatrix');

        examples_dir = fullfile(project_root, ...
            'geopdes_hierarchical_plasticity', 'examples', 'plasticity');
        old_dir = pwd;
        cleanup = onCleanup(@() cd(old_dir));
        cd(examples_dir);

        % Problem setup — small mesh, low refinement for speed
        problem_data.geo_name = 'geo_eighth_sphere.txt';
        problem_data.nmnn_sides   = [];
        problem_data.drchlt_sides = [];
        problem_data.press_sides  = [1];
        problem_data.symm_sides   = [3 4 5];
        problem_data.slider_sides = [6];
        problem_data.penalty_slider = @(x,y,z) 1e8*ones(size(x));
        E  = 210000; nu = 0.3;
        problem_data.yield_stress = @(x,y,z) 240*ones(size(x));
        problem_data.kappa_lame   = @(x,y,z) E/(3*(1-2*nu))*ones(size(x));
        problem_data.mu_lame      = @(x,y,z) E/(2*(1+nu))*ones(size(x));
        nload = 2;
        P = 332*0.99;
        problem_data.f = @(x,y,z) zeros(3,size(x,1),size(x,2),size(x,3));
        problem_data.g = @(x,y,z,ind) zeros(3,size(x,1),size(x,2),size(x,3));
        problem_data.h = @(x,y,z,ind) zeros(3,size(x,1),size(x,2),size(x,3));
        problem_data.p = @(x,y,z) P*ones(size(x));

        p = 2;
        method_data.degree         = [p p p];
        method_data.regularity     = [p-1 p-1 p-1];
        method_data.nsub_coarse    = [3, 1, 1];
        method_data.nsub_refine    = [2 2 2];
        method_data.nquad          = [p+1 p+1 p+1];
        method_data.space_type     = 'standard';
        method_data.truncated      = 1;
        method_data.nload          = nload;
        method_data.newton_tol     = 1e-6;
        method_data.newton_tol_abs = 1e-8;
        method_data.newton_iter_max = 50;
        method_data.type_projection = 'QI';

        adaptivity_data.flag        = 'elements';
        adaptivity_data.estimator   = 'sphere_front';
        adaptivity_data.C0_est      = 1.0;
        adaptivity_data.mark_param  = 0.5;
        adaptivity_data.mark_param_coarsening = 0.0;
        adaptivity_data.mark_strategy = 'MS';
        adaptivity_data.max_level   = 2;
        adaptivity_data.max_ndof    = 10000;
        adaptivity_data.num_max_iter = 2;
        adaptivity_data.max_nel     = 5000;
        adaptivity_data.tol         = 1e-5 * 3.14 * 30000;
        adaptivity_data.adm_strategy = 'admissible';
        adaptivity_data.coarsening_flag = 'any';
        adaptivity_data.adm         = p;

        % Run the solver
        [geometry, cell_hmsh, cell_hspace, ~, cell_u, ~, ~, ~] = ...
            adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);

        % Sanity checks
        ndof = cell_hspace{nload}.ndof;
        nel  = cell_hmsh{nload}.nel;
        assert(ndof > 0, 'Final ndof should be positive, got %d', ndof);
        assert(nel  > 0, 'Final nel should be positive, got %d', nel);
        assert(~any(isnan(cell_u{nload})), 'Solution contains NaN');
        assert(~any(isinf(cell_u{nload})), 'Solution contains Inf');

        R.passed  = true;
        R.message = sprintf('Completed: ndof=%d, nel=%d, 2 load steps', ndof, nel);
    catch ME
        R.passed  = false;
        R.message = ME.message;
    end
    R.time = toc(t0);
end

% ========================================================================
%  Suite: Benchmark (L2 vs QI vs QI_C0 accuracy comparison)
% ========================================================================
function R = run_benchmark(project_root)
    R.name = 'Benchmark: L2 vs QI vs QI_C0 (4x4x4, 3 load steps)';
    t0 = tic;
    try
        warning('off', 'MATLAB:nearlySingularMatrix');
        warning('off', 'MATLAB:singularMatrix');

        examples_dir = fullfile(project_root, ...
            'geopdes_hierarchical_plasticity', 'examples', 'plasticity');
        old_dir = pwd;
        cleanup = onCleanup(@() cd(old_dir));
        cd(examples_dir);

        % Run the benchmark script
        benchmark_qi_optimized;

        R.passed  = true;
        R.message = 'Benchmark completed successfully';
    catch ME
        R.passed  = false;
        R.message = ME.message;
    end
    R.time = toc(t0);
end
