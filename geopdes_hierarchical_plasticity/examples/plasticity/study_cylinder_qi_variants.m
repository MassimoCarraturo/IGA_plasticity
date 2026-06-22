% STUDY_CYLINDER_QI_VARIANTS  Convergence study for projection
% variants on the pressurised plane-strain thick-walled cylinder
% (quarter annulus, Hill's benchmark).  2D analogue of
% study_sphere_qi_variants.m.
%
% Six projection variants, each run for max_level = 1,2,3,4:
%   QI_C0      — standard (p=2) QI space with C^0 regularity
%   QI_p3      — cubic (p=3) QI space, C^2 regularity
%   QI_graded  — graded continuity: coarse levels C^1, fine levels C^0
%   L2         — L^2 projection onto the scalar stress space
%   QI         — plain QI projection (p=2, C^1 regularity)
%   QI_fine    — plain QI projection on a once-bisected scalar mesh
%                (num_bisections = 1) — feasible in 2D
%
% Uses the direct div(sigma) residual estimator (no C1 re-projection).
% Errors are compared against the Hill analytical plane-strain solution.
%
% Output:
%   results/cylinder_qi_study/  — .dat tables and standalone pgfplots .tex

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
results_dir = fullfile(here, 'results', 'cylinder_qi_study');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

addpath(genpath(fullfile(project_root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes-3.2.2', 'geopdes', 'inst')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'initialize')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
addpath(genpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
addpath(fullfile(project_root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

cd(here);
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

%% ========================================================================
%  PHYSICAL DATA
%  ========================================================================
E  = 210000;  nu = 0.3;  sigma_y = 240;        % [MPa]
a  = 100;     b  = 200;                         % inner / outer radius [mm]
% Pressure placing the elastic-plastic front at r = 160 mm (partially
% plastic, elastic outer zone retained).  Cylinder P_lim = 2*sy/sqrt(3)*ln(b/a).
pos_plastic_front = 160;
P  = 2*sigma_y/sqrt(3) * ( log(pos_plastic_front/a) + 0.5*(1 - pos_plastic_front^2/b^2) );
nload = 5;

problem_data.geo_name = 'geo_ring_SouzaNeto.txt';
problem_data.nmnn_sides   = [1];       % inner arc: pressure (Neumann traction g)
problem_data.drchlt_sides = [2];       % outer arc: exact-displacement Dirichlet
problem_data.press_sides  = [];
problem_data.symm_sides   = [3 4];     % two radial symmetry edges
problem_data.slider_sides = [];

problem_data.yield_stress = @(x, y) sigma_y * ones(size(x));
problem_data.kappa_lame   = @(x, y) E/(3*(1-2*nu)) * ones(size(x));
problem_data.mu_lame      = @(x, y) E/(2*(1+nu)) * ones(size(x));

problem_data.f = @(x, y) zeros(2, size(x,1), size(x,2));
problem_data.g = @(x, y, ind) test_plane_strain_ring_plasticity_g_nmnn(x, y, P, nu, ind);
problem_data.h = @(x, y, ind, mult) test_plane_strain_ring_plasticity_uex(x, y, E, nu, sigma_y, P, mult, ind);
problem_data.p = @(x, y, ind) P * ones(size(x));

problem_data.R_i = a;
problem_data.R_o = b;
problem_data.s_y = sigma_y;
problem_data.Pmax = P;

%% Hill analytical reference (plane strain)
n_points_ref = 200;
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius_ex] = ...
    Hill_solution(E, nu, sigma_y, a, b, P, nload, n_points_ref);
sigma_ex_fun = cell(1, nload);
P_steps = linspace(P/nload, P, nload);
for iload = 1:nload
    P_i = P_steps(iload);
    sigma_ex_fun{iload} = @(x, y) eval_Hill_cylinder(sigma_y, x, y, a, b, P_i);
end

%% ========================================================================
%  DISCRETIZATION TEMPLATE
%  ========================================================================
p = 2;
md0.degree         = [p p];
md0.regularity     = [p-1 p-1];
md0.nsub_coarse    = [2 2];
md0.nsub_refine    = [2 2];
md0.nquad          = [p+1 p+1];
md0.space_type     = 'standard';
md0.truncated      = 1;
md0.nload          = nload;
md0.newton_tol     = 1e-8;
md0.newton_tol_abs = 1e-10;
md0.newton_iter_max = 100;

ad0.flag                   = 'elements';
ad0.estimator              = 'div_sigma_direct';
ad0.C0_est                 = 1.0;
ad0.mark_param             = 0.5;
ad0.mark_param_coarsening  = 0.1;
ad0.mark_strategy          = 'MS';
ad0.max_ndof               = 15000;
ad0.max_nel                = 5000;
ad0.num_max_iter           = 3;
ad0.tol                    = 1e-5;   % scalar: max_level drives convergence, avoids premature stop
ad0.adm_strategy           = 'admissible';
ad0.coarsening_flag        = 'any';
ad0.adm                    = p;

%% ========================================================================
%  DEFINE CONFIGURATIONS: 6 methods × 4 adaptivity levels
%  ========================================================================
method_tags  = {'QI_C0', 'QI_p3', 'QI_graded', 'L2', 'QI', 'QI_fine'};
max_levels   = [1 2 3 4];
n_methods    = numel(method_tags);
n_levels     = numel(max_levels);

load_step_eval = [2 5];      % load steps at which to evaluate stress profiles

% Pre-allocate results (or load existing partial results)
checkpoint_file = fullfile(results_dir, 'cylinder_qi_checkpoint.mat');
if exist(checkpoint_file, 'file')
    tmp = load(checkpoint_file, 'R');
    R = tmp.R;
    fprintf('  Loaded checkpoint with existing results.\n');
else
    R = struct();
end
for im = 1:n_methods
    for il = 1:n_levels
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
        if ~isfield(R, tag) || R.(tag).ndof == 0
            R.(tag).ndof       = 0;
            R.(tag).nel        = 0;
            R.(tag).nlevels    = 0;
            R.(tag).wall       = 0;
            R.(tag).ndof_scalar = 0;
            R.(tag).err_sigma_r = nan(1, nload);
            R.(tag).err_sigma_t = nan(1, nload);
            R.(tag).err_sigma_r_rel = nan(1, nload);
            R.(tag).err_sigma_t_rel = nan(1, nload);
            R.(tag).u_r = zeros(1, nload+1);
            R.(tag).P_i = zeros(1, nload+1);
            R.(tag).sigma_rad  = [];
            R.(tag).sigma_tan  = [];
        end
    end
end

%% ========================================================================
%  RUN ALL CONFIGURATIONS
%  ========================================================================
for im = 1:n_methods
    for il = 1:n_levels
        method_data      = md0;
        adaptivity_data  = ad0;
        adaptivity_data.max_level = max_levels(il);

        switch method_tags{im}
            case 'QI_p3'
                method_data.type_projection = 'QI';
                method_data.degree_projection = [3 3];
                method_data.nquad = [4 4];     % sufficient for cubic basis

            case 'QI_C0'
                method_data.type_projection = 'QI_C0';

            case 'QI_graded'
                method_data.type_projection = 'QI_graded';
                % Graded continuity: levels 1-2 keep C^{p-1}, levels 3+ use C^0.
                method_data.graded_transition_level = 2;

            case 'L2'
                method_data.type_projection = 'L2';

            case 'QI'
                method_data.type_projection = 'QI';
                % Plain QI: same degree as displacement (p=2), C^1 regularity

            case 'QI_fine'
                method_data.type_projection = 'QI';
                method_data.num_bisections = 1;   % once-bisected scalar mesh
        end

        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));

        % Skip already-completed configurations
        if R.(tag).ndof > 0
            fprintf('\n  [SKIP] %s (max_level=%d) — already completed.\n', ...
                    method_tags{im}, max_levels(il));
            continue
        end

        fprintf('\n================================================================\n');
        fprintf('  %s  (max_level = %d)\n', method_tags{im}, max_levels(il));
        fprintf('================================================================\n');

        wall0 = tic;
        [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
         cell_u, ~, cell_sigma, solution_data, cell_hmsh_scalar] = ...
            adaptivity_J2_plasticity(problem_data, method_data, adaptivity_data);
        wall_total = toc(wall0);

        % ---- Summary numbers ------------------------------------------------
        hmsh_end   = cell_hmsh{end};
        hspace_end = cell_hspace{end};
        hspace_s   = cell_hspace_scalar{end};

        R.(tag).ndof       = hspace_end.ndof;
        R.(tag).nel        = hmsh_end.nel;
        R.(tag).nlevels    = hmsh_end.nlevels;
        R.(tag).wall       = wall_total;
        R.(tag).ndof_scalar = hspace_s.ndof;

        fprintf('  wall=%.1fs  ndof=%d  nel=%d  ndof_scalar=%d  nlevels=%d\n', ...
                wall_total, hspace_end.ndof, hmsh_end.nel, hspace_s.ndof, hmsh_end.nlevels);

        % ---- Load-displacement curve -----------------------------------------
        P_i = [0, P_steps];
        u_r = zeros(1, nload+1);
        for i = 1:nload
            [eu, ~] = sp_eval(cell_u{i}, cell_hspace{i}, geometry, {1, .5});
            u_r(i+1) = norm(eu);
        end
        R.(tag).u_r = u_r;
        R.(tag).P_i = P_i;

        % ---- L2 error of sigma_r / sigma_t vs Hill --------------------------
        for iload = 1:nload
            [errl2_r, errl2_t] = hspace_l2_error_stress( ...
                cell_hspace_scalar{iload}, cell_hmsh_scalar{iload}, ...
                cell_sigma{iload}, sigma_ex_fun{iload});
            [errl2_r_ex, errl2_t_ex] = hspace_l2_error_stress( ...
                cell_hspace_scalar{iload}, cell_hmsh_scalar{iload}, ...
                0*cell_sigma{iload}, sigma_ex_fun{iload});

            R.(tag).err_sigma_r(iload) = errl2_r;
            R.(tag).err_sigma_t(iload) = errl2_t;
            R.(tag).err_sigma_r_rel(iload) = errl2_r / errl2_r_ex;
            R.(tag).err_sigma_t_rel(iload) = errl2_t / errl2_t_ex;
        end

        % ---- Stress profiles along radial line (theta = 45 deg, eta = 0.5) ---
        voigt_inplane = [1, 2, 4];           % xx, yy, xy
        rad_pos = linspace(0, 1, length(radius_ex));

        sigma_rad = zeros(length(radius_ex), length(load_step_eval));
        sigma_tan = zeros(length(radius_ex), length(load_step_eval));
        for jl = 1:length(load_step_eval)
            lstep = load_step_eval(jl);
            for jp = 1:length(rad_pos)
                position = {rad_pos(jp), 0.5};
                s = zeros(1, 3);
                for k = 1:3
                    ic = voigt_inplane(k);
                    s(k) = sp_eval(cell_sigma{lstep}(:,ic), ...
                                   cell_hspace_scalar{lstep}, geometry, position);
                end
                % theta = 45 deg ray: sigma_r = (sxx+syy)/2 + sxy,
                %                     sigma_t = (sxx+syy)/2 - sxy
                sigma_rad(jp, jl) = 0.5*(s(1)+s(2)) + s(3);
                sigma_tan(jp, jl) = 0.5*(s(1)+s(2)) - s(3);
            end
        end
        R.(tag).sigma_rad = sigma_rad;
        R.(tag).sigma_tan = sigma_tan;

        fprintf('  err_sigma_r (final): abs=%.4e  rel=%.4e\n', ...
                R.(tag).err_sigma_r(end), R.(tag).err_sigma_r_rel(end));
        fprintf('  err_sigma_t (final): abs=%.4e  rel=%.4e\n', ...
                R.(tag).err_sigma_t(end), R.(tag).err_sigma_t_rel(end));

        % Incremental checkpoint save
        save(checkpoint_file, 'R');
        fprintf('  [checkpoint saved]\n');
    end
end

%% ========================================================================
%  SAVE RAW RESULTS
%  ========================================================================
save(fullfile(results_dir, 'cylinder_qi_results.mat'), 'R', ...
     'method_tags', 'max_levels', 'load_step_eval', ...
     'radius_ex', 'sigma_r_ex', 'sigma_t_ex', 'u_ex', 'P_ex', 'P_steps');

%% ========================================================================
%  CONSOLE SUMMARY TABLE
%  ========================================================================
fprintf('\n\n');
fprintf('========================================================================\n');
fprintf('  CYLINDER QI VARIANT STUDY — SUMMARY (final load step)\n');
fprintf('========================================================================\n');
fprintf('  %-12s  %4s  %6s  %6s  %6s  %10s  %10s  %8s\n', ...
        'Config', 'Lev', 'DOFs', 'Nel', 'DOFs_s', 'err_sig_r', 'err_sig_t', 'Wall[s]');
fprintf('------------------------------------------------------------------------\n');
for im = 1:n_methods
    for il = 1:n_levels
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
%  EXPORT DATA TABLES FOR pgfplots
%  ========================================================================

% --- Table 1: convergence data (error vs ndof for each method) -----------
for im = 1:n_methods
    fname = fullfile(results_dir, sprintf('convergence_%s.dat', method_tags{im}));
    fid = fopen(fname, 'w');
    fprintf(fid, 'max_level\tndof\tndof_scalar\tnel\terr_sigma_r\terr_sigma_t\terr_sigma_r_rel\terr_sigma_t_rel\twall\n');
    for il = 1:n_levels
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(il));
        fprintf(fid, '%d\t%d\t%d\t%d\t%.10e\t%.10e\t%.10e\t%.10e\t%.4f\n', ...
            max_levels(il), R.(tag).ndof, R.(tag).ndof_scalar, R.(tag).nel, ...
            R.(tag).err_sigma_r(end), R.(tag).err_sigma_t(end), ...
            R.(tag).err_sigma_r_rel(end), R.(tag).err_sigma_t_rel(end), ...
            R.(tag).wall);
    end
    fclose(fid);
end

% --- Table 2: load-displacement curves -----------------------------------
write_dat(fullfile(results_dir, 'load_disp_hill.dat'), {'u', 'P'}, [u_ex(:), P_ex(:)]);
for im = 1:n_methods
    tag = sprintf('%s_lev%d', method_tags{im}, max_levels(end));
    write_dat(fullfile(results_dir, sprintf('load_disp_%s.dat', method_tags{im})), ...
              {'u', 'P'}, [R.(tag).u_r(:), R.(tag).P_i(:)]);
end

% --- Table 3: stress profiles at selected load steps ---------------------
for jl = 1:length(load_step_eval)
    step = load_step_eval(jl);

    % sigma_r
    headers = {'r', 'hill'};
    cols = [radius_ex(:), sigma_r_ex(step,:).'];
    for im = 1:n_methods
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(end));
        headers{end+1} = method_tags{im}; %#ok<SAGROW>
        cols = [cols, R.(tag).sigma_rad(:, jl)]; %#ok<AGROW>
    end
    write_dat(fullfile(results_dir, sprintf('sigma_r_step%d.dat', step)), headers, cols);

    % sigma_t
    headers = {'r', 'hill'};
    cols = [radius_ex(:), sigma_t_ex(step,:).'];
    for im = 1:n_methods
        tag = sprintf('%s_lev%d', method_tags{im}, max_levels(end));
        headers{end+1} = method_tags{im}; %#ok<SAGROW>
        cols = [cols, R.(tag).sigma_tan(:, jl)]; %#ok<AGROW>
    end
    write_dat(fullfile(results_dir, sprintf('sigma_t_step%d.dat', step)), headers, cols);
end

% --- Table 4: per-load-step error evolution (for max_level=4) ------------
for im = 1:n_methods
    tag = sprintf('%s_lev%d', method_tags{im}, max_levels(end));
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
write_pgfplots(fullfile(results_dir, 'cylinder_qi_plots.tex'), ...
               method_tags, max_levels, load_step_eval, results_dir);

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
end


function [errl2_rad, errl2_tan] = hspace_l2_error_stress(hspace, hmsh, sigma, sigma_ex)
    errl2_rad = 0;  errl2_tan = 0;
    last_dof = cumsum(hspace.ndof_per_level);
    for ilev = 1:hmsh.nlevels
        if (hmsh.nel_per_level(ilev) > 0)
            msh_lev = msh_evaluate_element_list(hmsh.mesh_of_level(ilev), hmsh.active{ilev});
            sp_level = sp_evaluate_element_list(hspace.space_of_level(ilev), ...
                msh_lev, 'value', true);
            [e_r, e_t] = sp_l2_error_stress_2d(sp_level, msh_lev, ...
                hspace.Csub{ilev}*sigma(1:last_dof(ilev),:), sigma_ex);
            errl2_rad = errl2_rad + e_r^2;
            errl2_tan = errl2_tan + e_t^2;
        end
    end
    errl2_rad = sqrt(errl2_rad);
    errl2_tan = sqrt(errl2_tan);
end


function [errl2_rad, errl2_tan] = sp_l2_error_stress_2d(sp, msh, sigma, sigma_ex)
    % In-plane Voigt components: 1 = xx, 2 = yy, 4 = xy.
    x = cell(msh.rdim, 1);
    for idim = 1:msh.rdim
        x{idim} = reshape(msh.geo_map(idim,:,:), msh.nqn*msh.nel, 1);
    end
    [sigma_rad_ex, sigma_tan_ex] = feval(sigma_ex, x{1}, x{2});
    sigma_rad_ex = reshape(sigma_rad_ex, 1, msh.nqn*msh.nel);
    sigma_tan_ex = reshape(sigma_tan_ex, 1, msh.nqn*msh.nel);
    w = msh.quad_weights .* msh.jacdet;

    sxx = reshape(sp_eval_msh(sigma(:,1), sp, msh), 1, msh.nqn*msh.nel);
    syy = reshape(sp_eval_msh(sigma(:,2), sp, msh), 1, msh.nqn*msh.nel);
    sxy = reshape(sp_eval_msh(sigma(:,4), sp, msh), 1, msh.nqn*msh.nel);

    th = atan2(x{2}, x{1}).';
    c = cos(th);  s = sin(th);
    sigma_rad = sxx.*c.^2 + syy.*s.^2 + 2*sxy.*s.*c;
    sigma_tan = sxx.*s.^2 + syy.*c.^2 - 2*sxy.*s.*c;

    sigma_rad = reshape(sigma_rad, [1, msh.nqn, msh.nel]);
    sigma_tan = reshape(sigma_tan, [1, msh.nqn, msh.nel]);
    sigma_rad_ex = reshape(sigma_rad_ex, [1, msh.nqn, msh.nel]);
    sigma_tan_ex = reshape(sigma_tan_ex, [1, msh.nqn, msh.nel]);

    errl2_rad = sqrt(sum(sum(reshape(sum((sigma_rad-sigma_rad_ex).^2,1), [msh.nqn,msh.nel]) .* w)));
    errl2_tan = sqrt(sum(sum(reshape(sum((sigma_tan-sigma_tan_ex).^2,1), [msh.nqn,msh.nel]) .* w)));
end


function [sigma_r, sigma_t] = eval_Hill_cylinder(s_y, x, y, a, b, P)
    % Plane-strain von Mises thick-walled cylinder (Hill).  Y = 2*s_y/sqrt(3).
    Y   = 2*s_y/sqrt(3);
    P_0 = Y/2 * (1 - a^2/b^2);
    radius = sqrt(x.^2 + y.^2);
    sigma_r = zeros(size(radius));
    sigma_t = zeros(size(radius));
    if P < P_0
        % Fully elastic (Lamé)
        for j = 1:size(radius,2)
            for i = 1:size(radius,1)
                r = radius(i,j);
                sigma_r(i,j) = (P*a^2)/(b^2-a^2) - a^2*b^2*P/((b^2-a^2)*r^2);
                sigma_t(i,j) = (P*a^2)/(b^2-a^2) + a^2*b^2*P/((b^2-a^2)*r^2);
            end
        end
    else
        % Elastic-plastic: front radius c
        fun_front = @(xx) -P/Y + log(xx/a) + 0.5*(1 - xx.^2/b^2);
        c = fsolve(fun_front, 0.5*(a+b), optimoptions('fsolve','Display','off'));
        if ~(a < c && c < b)
            error('eval_Hill_cylinder: front radius c=%.3f outside (%.0f, %.0f) for P/Y=%.3f', c, a, b, P/Y);
        end
        for j = 1:size(radius,2)
            for i = 1:size(radius,1)
                r = radius(i,j);
                if r < c
                    sigma_r(i,j) = Y*(-0.5 - log(c/r) + c^2/(2*b^2));
                    sigma_t(i,j) = Y*( 0.5 - log(c/r) + c^2/(2*b^2));
                else
                    sigma_r(i,j) = -Y * c^2/(2*b^2) * (b^2/r^2 - 1);
                    sigma_t(i,j) =  Y * c^2/(2*b^2) * (b^2/r^2 + 1);
                end
            end
        end
    end
end


function write_pgfplots(filename, method_tags, max_levels, load_steps, results_dir)
    n_methods = numel(method_tags);

    % pgfplots style per method (must match method_tags order)
    styles = {
        'green!50!black,  mark=diamond*,  mark size=2pt, very thick'
        'blue!70!black,   mark=*,         mark size=2pt, very thick'
        'red!70!black,    mark=triangle*, mark size=2pt, very thick'
        'orange!80!black, mark=square*,   mark size=2pt, very thick'
        'violet,          mark=pentagon*, mark size=2pt, very thick'
        'cyan!60!black,   mark=oplus*,    mark size=2pt, very thick'};
    labels = {'QI ($C^0$)', 'QI ($p{=}3$, $C^2$)', 'QI (graded)', ...
              '$L^2$ proj.', 'QI (plain)', 'QI (finer mesh)'};

    fid = fopen(filename, 'w');
    fp = @(varargin) fprintf(fid, varargin{:});

    fp('%% Auto-generated by study_cylinder_qi_variants.m\n');
    fp('%% Compile: pdflatex cylinder_qi_plots.tex\n');
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
    fp('%%%% Load--displacement curves (max_level = %d)\n', max_levels(end));
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
        fp('\\addlegendentry{%s}\n', labels{im});
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
            fp('\\addplot[%s, mark repeat=8] table[x=r, y=%s] {sigma_r_step%d.dat};\n', ...
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
            fp('\\addplot[%s, mark repeat=8] table[x=r, y=%s] {sigma_t_step%d.dat};\n', ...
               styles{im}, method_tags{im}, step);
            fp('\\addlegendentry{%s}\n', labels{im});
        end
        fp('\\end{axis}\n');
        fp('\\end{tikzpicture}\n\n');
    end

    % ==== Plot 6: Error vs load step (max_level=4) =========================
    fp('%%%% Error evolution over load steps (max_level = %d)\n', max_levels(end));
    fp('\\begin{tikzpicture}\n');
    fp('\\begin{groupplot}[\n');
    fp('  group style={group size=2 by 1, horizontal sep=2cm},\n');
    fp('  width=8cm, height=7cm,\n');
    fp('  xlabel={Load step},\n');
    fp('  ymode=log,\n');
    fp('  legend pos=north west]\n');

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
end
