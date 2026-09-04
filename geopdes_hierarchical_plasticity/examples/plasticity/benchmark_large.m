% BENCHMARK_LARGE  Large-scale old-vs-new QI + L2 comparison.
%
%   Runs the eighth-sphere J2-plasticity benchmark on a 10x4x4 mesh with
%   up to 5 adaptive levels.  Compares:
%     1) L2 projection            (baseline — no libqi)
%     2) QI  with old libqi v1.0  (brute-force O(N) scan)
%     3) QI  with new libqi v1.1  (spatial index + all optimisations)
%     4) QI_C0 with old libqi v1.0
%     5) QI_C0 with new libqi v1.1
%
%   Swaps MEX binaries between old/new builds.  Old binaries must exist as
%   qi_mex_old.mexw64 / qi_local_ls_mex_old.mexw64 alongside the current ones.

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
% Problem setup — LARGE mesh
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

p = 2;
method_data_template.degree      = [p p p];
method_data_template.regularity  = [p-1 p-1 p-1];
method_data_template.nsub_coarse = [10, 4, 4];        % <-- LARGER initial mesh
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
adaptivity_data.max_level   = 5;                      % <-- MORE levels
adaptivity_data.max_ndof    = 30000;                   % <-- LARGER cap
adaptivity_data.num_max_iter = 3;                      % <-- MORE refinement iters
adaptivity_data.max_nel     = 15000;
adaptivity_data.tol         = 1e-5 * 3.14 * 30000;
adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any';
adaptivity_data.adm         = p;

% ========================================================================
% Define runs
% ========================================================================
runs = {
%   label         method   mex_version
    'L2',         'L2',    'new'   % L2 doesn't use libqi, version is irrelevant
    'QI_old',     'QI',    'old'
    'QI_new',     'QI',    'new'
    'QI_C0_old',  'QI_C0', 'old'
    'QI_C0_new',  'QI_C0', 'new'
};

results = struct();

for irun = 1:size(runs, 1)
    label       = runs{irun, 1};
    method_name = runs{irun, 2};
    mex_ver     = runs{irun, 3};

    fprintf('\n================================================================\n');
    fprintf('  %s : type_projection = %s, libqi = %s\n', label, method_name, mex_ver);
    fprintf('================================================================\n');

    % Swap MEX if needed
    swap_mex(mex_dir, mex_ver);

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

    % L2 stress error at last load step
    P_i_last = P / nload * nload;   % = P
    sigma_ex = @(x,y,z) eval_Hill_sphere(s_y, x,y,z, 100, 200, P_i_last);
    [errl2_rad, errl2_tan] = hspace_l2_error_stress( ...
        cell_hspace_scalar{nload}, cell_hmsh{nload}, cell_sigma{nload}, sigma_ex);
    [errl2_rad_ex, errl2_tan_ex] = hspace_l2_error_stress( ...
        cell_hspace_scalar{nload}, cell_hmsh{nload}, 0*cell_sigma{nload}, sigma_ex);

    R = struct();
    R.label      = label;
    R.method     = method_name;
    R.mex_ver    = mex_ver;
    R.wall_time  = t_wall;
    R.ndof       = ndof_final;
    R.nel        = nel_final;
    R.errl2_rad  = errl2_rad;
    R.errl2_tan  = errl2_tan;
    R.errl2_rad_rel = errl2_rad / errl2_rad_ex;
    R.errl2_tan_rel = errl2_tan / errl2_tan_ex;
    results.(matlab.lang.makeValidName(label)) = R;

    fprintf('\n  %s done: wall = %.2f s, ndof = %d, nel = %d\n', ...
            label, t_wall, ndof_final, nel_final);
    fprintf('  L2 err sigma_r: abs=%.4e  rel=%.4e\n', errl2_rad, errl2_rad/errl2_rad_ex);
    fprintf('  L2 err sigma_t: abs=%.4e  rel=%.4e\n', errl2_tan, errl2_tan/errl2_tan_ex);
end

% Restore new MEX at exit
swap_mex(mex_dir, 'new');

% ========================================================================
% Summary table
% ========================================================================
fprintf('\n\n');
fprintf('================================================================================\n');
fprintf('  LARGE-SCALE BENCHMARK — libqi v1.0 (old) vs v1.1 (new) + L2\n');
fprintf('================================================================================\n');
fprintf('%-12s  %10s  %7s  %7s  %12s  %12s  %12s  %12s\n', ...
    'Run', 'WallTime', 'ndof', 'nel', 'L2_sr_abs', 'L2_sr_rel', 'L2_st_abs', 'L2_st_rel');
fprintf('%-12s  %10s  %7s  %7s  %12s  %12s  %12s  %12s\n', ...
    '---', '--------', '----', '---', '---------', '---------', '---------', '---------');

for irun = 1:size(runs, 1)
    label = runs{irun, 1};
    fn = matlab.lang.makeValidName(label);
    R = results.(fn);
    fprintf('%-12s  %9.2f s  %7d  %7d  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
        R.label, R.wall_time, R.ndof, R.nel, ...
        R.errl2_rad, R.errl2_rad_rel, R.errl2_tan, R.errl2_tan_rel);
end

% Speedup summary
fprintf('\nSpeedups:\n');
fn_qi_old  = matlab.lang.makeValidName('QI_old');
fn_qi_new  = matlab.lang.makeValidName('QI_new');
fn_c0_old  = matlab.lang.makeValidName('QI_C0_old');
fn_c0_new  = matlab.lang.makeValidName('QI_C0_new');
fn_l2      = matlab.lang.makeValidName('L2');

t_qi_old = results.(fn_qi_old).wall_time;
t_qi_new = results.(fn_qi_new).wall_time;
t_c0_old = results.(fn_c0_old).wall_time;
t_c0_new = results.(fn_c0_new).wall_time;
t_l2     = results.(fn_l2).wall_time;

fprintf('  QI  v1.1 vs v1.0 :  %.2f s -> %.2f s  = %.2fx speedup\n', t_qi_old, t_qi_new, t_qi_old/t_qi_new);
fprintf('  QI_C0 v1.1 vs v1.0: %.2f s -> %.2f s  = %.2fx speedup\n', t_c0_old, t_c0_new, t_c0_old/t_c0_new);
fprintf('  QI_new vs L2      : %.2f s vs %.2f s  = %.2fx\n', t_qi_new, t_l2, t_l2/t_qi_new);
fprintf('  QI_C0_new vs L2   : %.2f s vs %.2f s  = %.2fx\n', t_c0_new, t_l2, t_l2/t_c0_new);
fprintf('\n');

save(fullfile(results_dir, 'benchmark_large.mat'), 'results', 'runs');
fprintf('Results saved to %s\n', fullfile(results_dir, 'benchmark_large.mat'));

% ========================================================================
% Local functions
% ========================================================================
function swap_mex(mex_dir, tag)
    qi_active = fullfile(mex_dir, 'qi_mex.mexw64');
    qi_old    = fullfile(mex_dir, 'qi_mex_old.mexw64');
    qi_save   = fullfile(mex_dir, 'qi_mex_new.mexw64');
    ls_active = fullfile(mex_dir, 'qi_local_ls_mex.mexw64');
    ls_old    = fullfile(mex_dir, 'qi_local_ls_mex_old.mexw64');
    ls_save   = fullfile(mex_dir, 'qi_local_ls_mex_new.mexw64');

    clear qi_mex qi_local_ls_mex;

    if strcmpi(tag, 'old')
        if ~isfile(qi_save)
            copyfile(qi_active, qi_save);
            copyfile(ls_active, ls_save);
        end
        copyfile(qi_old, qi_active);
        copyfile(ls_old, ls_active);
    else
        if isfile(qi_save)
            copyfile(qi_save, qi_active);
            copyfile(ls_save, ls_active);
            delete(qi_save);
            delete(ls_save);
        end
    end
end

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
