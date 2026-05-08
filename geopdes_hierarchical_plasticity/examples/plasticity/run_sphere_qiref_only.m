% RUN_SPHERE_QIREF_ONLY  Re-run only the QI_ref pass for the eighth-sphere
% benchmark (with reduced num_bisections inside adaptivity_J2_plasticity).
% L2 and QI results are loaded from previously-saved .mat files.

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

% Suppress the "Matrix is close to singular" warnings that QI_ref emits a lot of
warning ('off', 'MATLAB:nearlySingularMatrix');

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
P = 332*.99;
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3));
problem_data.p = @(x, y, z) P*ones (size (x));

p = 2;
method_data.degree     = [p p p];
method_data.regularity = [p-1 p-1 p-1];
method_data.nsub_coarse = [15, 5, 5];
method_data.nsub_refine = [2 2 2];
method_data.nquad      = [p+1 p+1 p+1];
method_data.space_type = 'standard';
method_data.truncated  = 1;
method_data.nload      = nload;
method_data.newton_tol = 1e-8;
method_data.newton_tol_abs = 1e-10;
method_data.newton_iter_max = 100;
method_data.type_projection = 'QI_ref';

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

[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution (E, nu, problem_data.yield_stress(1), 100, 200, P, nload);
load_step_eval_stress = [2, 5];

fprintf('Running QI_ref...\n');
tic;
[geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
 cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
    adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);
fprintf('QI_ref solver done in %.1f s\n', toc);

% Post-process
u_r = zeros(1, method_data.nload+1);
P_i = zeros(1, method_data.nload+1);
for i = 1:method_data.nload+1
    P_i(i) = P/method_data.nload * (i-1);
    if i > 1
        eu = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1, .5, .5});
        u_r(i) = norm(eu);
    end
end

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

R = struct();
R.method      = 'QI_ref';
R.u_r         = u_r;
R.P_i         = P_i;
R.radius      = radius;
R.sigma_rad   = sigma_rad;
R.sigma_tan   = sigma_tan;
R.load_steps  = load_step_eval_stress;
R.solution    = solution_data;
save (fullfile(here, 'sphere_results_QI_ref.mat'), '-struct', 'R');

fprintf('QI_REF DONE\n');
