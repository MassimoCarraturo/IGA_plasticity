% BENCHMARK_QI_OPTIMIZED  Comparative benchmark: L2 vs QI vs QI_C0.
%
%   Runs the eighth-sphere J2-plasticity benchmark with a moderate mesh and
%   measures wall-clock time and L2 stress error for each projection method.
%   Uses the optimised libqi v1.1 (spatial index, symmetric S, power tables,
%   pre-computed knots, etc.).
%
%   Output:
%     - Console summary table: method, total time, projection time, L2 errors
%     - benchmark_results.mat saved to the current directory

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

% ========================================================================
% Problem setup (identical to ex_plastic_sphere_hier but with a practical mesh)
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

s_y = problem_data.yield_stress(1);
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution(E, nu, s_y, 100, 200, P, nload);

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
% Run each projection method with timing
% ========================================================================
proj_methods = {'L2', 'QI', 'QI_C0'};
load_step_eval = nload;   % evaluate stress error at last load step

bench = struct();

for imethod = 1:numel(proj_methods)
    method_name = proj_methods{imethod};
    fprintf('\n================================================================\n');
    fprintf('  BENCHMARK: type_projection = %s\n', method_name);
    fprintf('================================================================\n');

    method_data = method_data_template;
    method_data.type_projection = method_name;

    % Clear the persistent localLS function handle to force re-detection
    clear history_variable_projection_hier;

    % Wall-clock timing (total simulation)
    t_wall_start = tic;

    [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
     cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
        adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);

    t_wall = toc(t_wall_start);

    % ---- Post-processing: load-displacement curve --------------------------
    u_r = zeros(1, nload+1);
    P_i = zeros(1, nload+1);
    for i = 1:nload+1
        P_i(i) = P/nload * (i-1);
        if i > 1
            [eu, ~] = sp_eval(cell_u{i-1}, cell_hspace{i-1}, geometry, {1, .5, .5});
            u_r(i) = norm(eu);
        else
            u_r(i) = 0;
        end
    end

    % ---- L2 stress error vs Hill's analytical solution ---------------------
    errl2_sigma_rad = 0;
    errl2_sigma_tan = 0;
    errl2_sigma_rad_ex = 0;
    errl2_sigma_tan_ex = 0;

    k = load_step_eval;
    sigma_ex_fun = @(x, y, z) eval_Hill_sphere(s_y, x, y, z, 100, 200, P_i(k + 1));

    [errl2_sigma_rad, errl2_sigma_tan] = hspace_l2_error_stress( ...
        cell_hspace_scalar{k}, cell_hmsh{k}, cell_sigma{k}, sigma_ex_fun);
    [errl2_sigma_rad_ex, errl2_sigma_tan_ex] = hspace_l2_error_stress( ...
        cell_hspace_scalar{k}, cell_hmsh{k}, 0*cell_sigma{k}, sigma_ex_fun);

    % ---- Mesh statistics ---------------------------------------------------
    ndof_final  = cell_hspace{nload}.ndof;
    nel_final   = sum(cell_hmsh{nload}.nel_per_level);

    % ---- Store results -----------------------------------------------------
    B = struct();
    B.method         = method_name;
    B.wall_time      = t_wall;
    B.u_r            = u_r;
    B.P_i            = P_i;
    B.errl2_rad      = errl2_sigma_rad;
    B.errl2_tan      = errl2_sigma_tan;
    B.errl2_rad_rel  = errl2_sigma_rad / errl2_sigma_rad_ex;
    B.errl2_tan_rel  = errl2_sigma_tan / errl2_sigma_tan_ex;
    B.ndof_final     = ndof_final;
    B.nel_final      = nel_final;
    B.solution_data  = solution_data;

    bench.(matlab.lang.makeValidName(method_name)) = B;

    fprintf('\n  %s complete: wall = %.2f s, ndof = %d, nel = %d\n', ...
            method_name, t_wall, ndof_final, nel_final);
    fprintf('  L2 error (sigma_r): abs = %.4e, rel = %.4e\n', ...
            errl2_sigma_rad, errl2_sigma_rad / errl2_sigma_rad_ex);
    fprintf('  L2 error (sigma_t): abs = %.4e, rel = %.4e\n', ...
            errl2_sigma_tan, errl2_sigma_tan / errl2_sigma_tan_ex);
end

% ========================================================================
% Summary table
% ========================================================================
fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  COMPARATIVE RESULTS (optimised libqi v1.1)\n');
fprintf('================================================================\n');
fprintf('%-8s  %10s  %8s  %8s  %12s  %12s  %12s  %12s\n', ...
        'Method', 'WallTime', 'ndof', 'nel', ...
        'L2_sr_abs', 'L2_sr_rel', 'L2_st_abs', 'L2_st_rel');
fprintf('%-8s  %10s  %8s  %8s  %12s  %12s  %12s  %12s\n', ...
        '------', '---------', '----', '---', ...
        '---------', '---------', '---------', '---------');

for imethod = 1:numel(proj_methods)
    mn = matlab.lang.makeValidName(proj_methods{imethod});
    B = bench.(mn);
    fprintf('%-8s  %9.2f s  %8d  %8d  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
            B.method, B.wall_time, B.ndof_final, B.nel_final, ...
            B.errl2_rad, B.errl2_rad_rel, B.errl2_tan, B.errl2_tan_rel);
end

% Speedup relative to L2
if isfield(bench, 'L2')
    fprintf('\nSpeedup relative to L2:\n');
    t_l2 = bench.L2.wall_time;
    for imethod = 1:numel(proj_methods)
        mn = matlab.lang.makeValidName(proj_methods{imethod});
        B = bench.(mn);
        fprintf('  %-8s: %.2fx\n', B.method, t_l2 / B.wall_time);
    end
end

fprintf('\n');

% ========================================================================
% Save
% ========================================================================
save(fullfile(results_dir, 'benchmark_results.mat'), 'bench', 'proj_methods', ...
     'u_ex', 'P_ex', 'sigma_r_ex', 'sigma_t_ex', 'radius');
fprintf('Results saved to %s\n', fullfile(results_dir, 'benchmark_results.mat'));

% ========================================================================
% Local functions (duplicated from ex_plastic_sphere_hier.m for standalone use)
% ========================================================================
function [errl2_rad, errl2_tan] = hspace_l2_error_stress(hspace, hmsh, sigma, sigma_ex)
    errl2_rad = 0; errl2_tan = 0;
    last_dof = cumsum(hspace.ndof_per_level);
    for ilev = 1:hmsh.nlevels
        if (hmsh.nel_per_level(ilev) > 0)
            msh_level = hmsh.msh_lev{ilev};
            sp_level = sp_evaluate_element_list(hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
            [errl2_rad_lev, errl2_tan_lev] = sp_l2_error_stress(sp_level, msh_level, ...
                hspace.Csub{ilev}*sigma(1:last_dof(ilev),:), sigma_ex);
            errl2_rad = errl2_rad + errl2_rad_lev.^2;
            errl2_tan = errl2_tan + errl2_tan_lev.^2;
        end
    end
    errl2_rad = sqrt(errl2_rad);
    errl2_tan = sqrt(errl2_tan);
end

function [errl2_rad, errl2_tan] = sp_l2_error_stress(sp, msh, sigma, sigma_ex)
    voigt = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];
    for idir = 1:msh.rdim
        x{idir} = reshape(msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
    end
    [sigma_rad_ex, sigma_tan_ex] = feval(sigma_ex, x{:});
    sigma_rad_ex = reshape(sigma_rad_ex, sp.ncomp, msh.nqn, msh.nel);
    sigma_tan_ex = reshape(sigma_tan_ex, sp.ncomp, msh.nqn, msh.nel);
    w = msh.quad_weights .* msh.jacdet;
    sigma_matrix = zeros(3,3,msh.nqn*msh.nel);
    for icomp = 1:6
        eu = sp_eval_msh(sigma(:,icomp), sp, msh);
        eu = reshape(eu, sp.ncomp, msh.nqn * msh.nel);
        sigma_matrix(voigt(icomp,1), voigt(icomp,2), :) = eu;
        if icomp > 3
            sigma_matrix(voigt(icomp,2), voigt(icomp,1), :) = eu;
        end
    end
    sigma_rad = zeros(1, msh.nqn*msh.nel);
    sigma_tan = zeros(1, msh.nqn*msh.nel);
    for iGP = 1:msh.nqn*msh.nel
        coordX = x{1}(iGP); coordY = x{2}(iGP); coordZ = x{3}(iGP);
        rad_dir = [coordX; coordY; coordZ] ./ norm([coordX; coordY; coordZ]);
        tan_dir = [coordY; -coordX; 0] ./ norm([coordX; coordY; 0]);
        sigma_rad(iGP) = rad_dir' * sigma_matrix(:,:,iGP) * rad_dir;
        sigma_tan(iGP) = tan_dir' * sigma_matrix(:,:,iGP) * tan_dir;
    end
    sigma_rad = reshape(sigma_rad, [1, msh.nqn, msh.nel]);
    sigma_tan = reshape(sigma_tan, [1, msh.nqn, msh.nel]);
    errl2_rad_elem = sum(reshape(sum((sigma_rad - sigma_rad_ex).^2, 1), [msh.nqn, msh.nel]) .* w);
    errl2_tan_elem = sum(reshape(sum((sigma_tan - sigma_tan_ex).^2, 1), [msh.nqn, msh.nel]) .* w);
    errl2_rad = sqrt(sum(errl2_rad_elem));
    errl2_tan = sqrt(sum(errl2_tan_elem));
end

function [sigma_r, sigma_t] = eval_Hill_sphere(s_y, x, y, z, a, b, P)
    P_0 = 2*s_y/3*(1-(a^3/b^3));
    radius = sqrt(x.^2 + y.^2 + z.^2);
    sigma_r = zeros(size(radius));
    sigma_t = zeros(size(radius));
    if P < P_0
        for j = 1:size(radius, 2)
            for i = 1:size(radius, 1)
                sigma_r(i,j) = -P*a^3 / (b^3 - a^3) * (b^3 / radius(i,j)^3 - 1);
                sigma_t(i,j) =  P*a^3 / (b^3 - a^3) * (0.5*b^3 / radius(i,j)^3 + 1);
            end
        end
    else
        fun_front = @(x) -P + 2*s_y*log(x/a) + 2/3*s_y*(1 - x.^3/b^3);
        c = fsolve(fun_front, 100, optimoptions('fsolve', 'Display', 'off'));
        for j = 1:size(radius, 2)
            for i = 1:size(radius, 1)
                if radius(i,j) <= c
                    sigma_r(i,j) = -2*s_y*(log(c/radius(i,j)) + 1/3*(1-c^3/b^3));
                    sigma_t(i,j) =  2*s_y*(0.5 - log(c/radius(i,j)) - 1/3*(1-c^3/b^3));
                else
                    sigma_r(i,j) = -2*s_y*c^3/(3*b^3) * (b^3/radius(i,j)^3 - 1);
                    sigma_t(i,j) =  2*s_y*c^3/(3*b^3) * (0.5*b^3/radius(i,j)^3 + 1);
                end
            end
        end
    end
end
