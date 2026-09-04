% EX_PLASTIC_SPHERE_HIER  Eighth-sphere J2-plasticity benchmark from
% "Computational Methods for Plasticity: Theory and Applications"
% (de Souza Neto, Peric, Owen).  Compares three projection methods for the
% history variables — L2, QI, QI_C0 — against Hill's analytical solution.
%
% Outputs written to results/:
%   results/sphere_results_<method>.mat     numerical curves for each method
%   results/pgfplots/sphere_load_disp.dat   load-displacement table (all methods + Hill)
%   results/pgfplots/sphere_sigma_r_<step>.dat   radial stress at chosen load step
%   results/pgfplots/sphere_sigma_t_<step>.dat   tangential stress at chosen load step
%   results/pgfplots/sphere_comparison.tex  standalone pgfplots document
%   results/sphere_sigma_r.png / sphere_sigma_t.png  MATLAB comparison figures

% 1) PHYSICAL DATA OF THE PROBLEM
clear all
clc
close all

here = fileparts(mfilename('fullpath'));
results_dir = fullfile(here, 'results');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_eighth_sphere.txt'; % direction 1=radial

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [1];
problem_data.symm_sides   = [3 4 5];
problem_data.slider_sides = [6]; % zero displacement
problem_data.penalty_slider = @(x, y, z) 1e8 * ones (size (x));

% Physical parameters
E  =  210000;                                  % MPa
nu = 0.3;                                      % -
sigma_y = 240;
problem_data.yield_stress = @(x, y, z) sigma_y * ones (size (x));  % MPa
problem_data.kappa_lame = @(x, y, z) E/(3*(1-2*nu)) * ones (size (x));
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
nload = 5;
P = 332*.99;                                        % Limit internal pressure [MPa]
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.p = @(x, y, z) P*ones (size (x));

% Hill (analytical) reference
s_y = problem_data.yield_stress(1);
[u_ex, P_ex, sigma_r_ex, sigma_t_ex, radius] = ...
    sphere_solution (E, nu, s_y, 100, 200, P, nload);
load_step_eval_stress = [2, 5];

% 2) DISCRETIZATION COMMON TO ALL METHODS
p = 2;
method_data_template.degree     = [p p p];
method_data_template.regularity = [p-1 p-1 p-1];
method_data_template.nsub_coarse = [2,2,2];

method_data_template.nsub_refine = [2 2 2];
method_data_template.nquad      = [p+1 p+1 p+1];
method_data_template.space_type  = 'standard';
method_data_template.truncated   = 1;
method_data_template.nload      = nload;
method_data_template.newton_tol = 1e-8;
method_data_template.newton_tol_abs = 1e-10;
method_data_template.newton_iter_max = 100;

adaptivity_data.flag = 'elements';
adaptivity_data.estimator = 'sphere_front';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = 0.9;
adaptivity_data.mark_param_coarsening = 0.1;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 4;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 4;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5 * 3.14 * 30000;
adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any';
adaptivity_data.adm = p;

% 3) RUN EACH PROJECTION METHOD
proj_methods = {'L2', 'QI', 'QI_C0'};
results = struct();

for imethod = 1:numel(proj_methods)
    method_name = proj_methods{imethod};
    fprintf('\n================================================================\n');
    fprintf('Running ex_plastic_sphere_hier with type_projection = %s\n', method_name);
    fprintf('================================================================\n');

    method_data = method_data_template;
    method_data.type_projection = method_name;

    [geometry, cell_hmsh, cell_hspace, cell_hspace_scalar, ...
     cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
        adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);

    % ---- Post-processing: load-displacement curve ---------------------------
    u_r = zeros(1, method_data.nload+1);
    P_i = zeros(1, method_data.nload+1);
    for i = 1:method_data.nload+1
        P_i(i) = P/method_data.nload * (i-1);
        if i > 1
            [eu, ~] = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1, .5, .5});
            u_r(i) = norm(eu);
        else
            u_r(i) = 0;
        end
    end

    % ---- Post-processing: radial / tangential stress profiles --------------
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

        sigma_ex = @(x, y, z) eval_Hill_sphere(s_y, x, y, z, 100, 200, P_i(load_step_eval_stress(jload) + 1));
        
        [errl2_sigma_rad(jload), errl2_sigma_tan(jload)] = hspace_l2_error_stress (cell_hspace_scalar{load_step_eval_stress(jload)}, ...
            cell_hmsh{load_step_eval_stress(jload)}, cell_sigma{load_step_eval_stress(jload)}, sigma_ex);
        [errl2_sigma_rad_ex(jload), errl2_sigma_tan_ex(jload)] = hspace_l2_error_stress (cell_hspace_scalar{load_step_eval_stress(jload)}, ...
            cell_hmsh{load_step_eval_stress(jload)}, 0*cell_sigma{load_step_eval_stress(jload)}, sigma_ex);
    end

    


    R = struct();
    R.method      = method_name;
    R.u_r         = u_r;
    R.P_i         = P_i;
    R.radius      = radius;
    R.sigma_rad   = sigma_rad;
    R.sigma_tan   = sigma_tan;
    R.errl2_sigma_rad   = errl2_sigma_rad;
    R.errl2_sigma_tan   = errl2_sigma_tan;
    R.errl2_rel_sigma_rad   = errl2_sigma_rad ./ errl2_sigma_rad_ex;
    R.errl2_rel_sigma_tan   = errl2_sigma_tan ./ errl2_sigma_tan_ex;
    R.load_steps  = load_step_eval_stress;
    R.solution    = solution_data;
    results.(matlab.lang.makeValidName(method_name)) = R;

    save (fullfile(results_dir, sprintf('sphere_results_%s.mat', method_name)), '-struct', 'R');
end

% 4) PLOT COMPARISON IN MATLAB
out_dir = fullfile(results_dir, 'pgfplots');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

method_styles = {'-o', '-s', '-d'};
method_colors = {[0.20 0.45 0.80], [0.85 0.40 0.20], [0.20 0.65 0.30]};
method_labels = proj_methods;

% Figure 1: load-displacement curve (Hill + three methods)
figure(1); clf; hold on; grid on; box on;
plot (u_ex, P_ex, '-k', 'LineWidth', 1.5, 'DisplayName', 'Hill (analytical)');
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    plot (R.u_r, R.P_i, method_styles{imethod}, ...
          'Color', method_colors{imethod}, 'LineWidth', 1.2, ...
          'MarkerSize', 6, 'MarkerFaceColor', method_colors{imethod}, ...
          'DisplayName', method_labels{imethod});
end
xlabel('u_r at outer boundary [mm]');
ylabel('Internal pressure P [MPa]');
title('Eighth-sphere: load–displacement curve');
legend('Location', 'northwest');

% Figure 2: radial stress profiles
figure(2); clf; hold on; grid on; box on;
for jload = 1:length(load_step_eval_stress)
    plot (radius, sigma_r_ex(load_step_eval_stress(jload),:), '-k', ...
          'LineWidth', 1.5, 'DisplayName', ...
          sprintf('Hill, step %d', load_step_eval_stress(jload)));
end
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    for jload = 1:length(load_step_eval_stress)
        plot (radius, R.sigma_rad(:, jload), method_styles{imethod}, ...
              'Color', method_colors{imethod}, ...
              'MarkerIndices', 1:10:length(radius), ...
              'DisplayName', sprintf('%s, step %d', ...
                                     method_labels{imethod}, ...
                                     load_step_eval_stress(jload)));
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_r [MPa]');
title('Eighth-sphere: radial stress');
legend('Location', 'best');

% Figure 3: tangential stress profiles
figure(3); clf; hold on; grid on; box on;
for jload = 1:length(load_step_eval_stress)
    plot (radius, sigma_t_ex(load_step_eval_stress(jload),:), '-k', ...
          'LineWidth', 1.5, 'DisplayName', ...
          sprintf('Hill, step %d', load_step_eval_stress(jload)));
end
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    for jload = 1:length(load_step_eval_stress)
        plot (radius, R.sigma_tan(:, jload), method_styles{imethod}, ...
              'Color', method_colors{imethod}, ...
              'MarkerIndices', 1:10:length(radius), ...
              'DisplayName', sprintf('%s, step %d', ...
                                     method_labels{imethod}, ...
                                     load_step_eval_stress(jload)));
    end
end
xlabel('Radial coordinate r [mm]');
ylabel('\sigma_t [MPa]');
title('Eighth-sphere: tangential stress');
legend('Location', 'best');

% Save MATLAB figures
saveas(figure(1), fullfile(results_dir, 'sphere_load_disp.png'));
saveas(figure(2), fullfile(results_dir, 'sphere_sigma_r.png'));
saveas(figure(3), fullfile(results_dir, 'sphere_sigma_t.png'));

% 5) EXPORT DATA TABLES FOR pgfplots
% Load–displacement table (one column per method + Hill)
ld_table = cell(0,1);
ld_header = sprintf('u_hill\tP_hill');
ld_rows_hill = [u_ex(:), P_ex(:)];
% pad / merge: each method shares the same P_i (load steps known a priori),
% so we write a separate file per series for simplicity.
write_dat (fullfile(out_dir, 'sphere_load_disp_hill.dat'), ...
           {'u', 'P'}, [u_ex(:), P_ex(:)]);
for imethod = 1:numel(proj_methods)
    R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
    write_dat (fullfile(out_dir, sprintf('sphere_load_disp_%s.dat', proj_methods{imethod})), ...
               {'u', 'P'}, [R.u_r(:), R.P_i(:)]);
end

% Stress profiles — one .dat per (quantity, load step), with a column per method + Hill
for jload = 1:length(load_step_eval_stress)
    step = load_step_eval_stress(jload);

    headers = {'r', 'hill'};
    cols    = [radius(:), sigma_r_ex(step,:).'];
    for imethod = 1:numel(proj_methods)
        R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
        headers{end+1} = proj_methods{imethod}; %#ok<SAGROW>
        cols = [cols, R.sigma_rad(:, jload)]; %#ok<AGROW>
    end
    write_dat (fullfile(out_dir, sprintf('sphere_sigma_r_step%d.dat', step)), headers, cols);

    headers = {'r', 'hill'};
    cols    = [radius(:), sigma_t_ex(step,:).'];
    for imethod = 1:numel(proj_methods)
        R = results.(matlab.lang.makeValidName(proj_methods{imethod}));
        headers{end+1} = proj_methods{imethod}; %#ok<SAGROW>
        cols = [cols, R.sigma_tan(:, jload)]; %#ok<AGROW>
    end
    write_dat (fullfile(out_dir, sprintf('sphere_sigma_t_step%d.dat', step)), headers, cols);
end

% 6) GENERATE pgfplots LaTeX SCRIPT
write_pgfplots_tex (fullfile(out_dir, 'sphere_comparison.tex'), ...
                    proj_methods, load_step_eval_stress);

fprintf('\nAll done. Data tables and LaTeX in %s/\n', out_dir);

% =========================================================================
% Local helper functions
% =========================================================================
function write_dat (filename, headers, data)
% Write a tab-separated ASCII table for pgfplots `table[header=true]`.
    fid = fopen (filename, 'w');
    if fid < 0; error('Cannot open %s for writing', filename); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf (fid, '%s', headers{1});
    for k = 2:numel(headers); fprintf (fid, '\t%s', headers{k}); end
    fprintf (fid, '\n');
    ncols = size(data, 2);
    fmt = [strjoin(repmat({'%.10g'}, 1, ncols), '\t') '\n'];
    for i = 1:size(data,1)
        fprintf (fid, fmt, data(i,:));
    end
end

function write_pgfplots_tex (filename, proj_methods, load_steps)
    method_pgf_styles = {
        'mark=*,        mark size=2pt, color=blue!70!black,    very thick'
        'mark=square*,  mark size=2pt, color=orange!80!black,  very thick'
        'mark=diamond*, mark size=2pt, color=green!50!black,   very thick'};
    fid = fopen (filename, 'w');
    if fid < 0; error('Cannot open %s for writing', filename); end
    cleanup = onCleanup(@() fclose(fid));

    fprintf (fid, '%% Auto-generated by ex_plastic_sphere_hier.m\n');
    fprintf (fid, '%% Compile with: pdflatex sphere_comparison.tex\n');
    fprintf (fid, '\\documentclass[tikz,border=4pt]{standalone}\n');
    fprintf (fid, '\\usepackage{pgfplots}\n');
    fprintf (fid, '\\pgfplotsset{compat=1.18, every axis/.append style={font=\\small}}\n\n');
    fprintf (fid, '\\begin{document}\n\n');

    % --- Load-displacement plot ----------------------------------------------
    fprintf (fid, '\\begin{tikzpicture}\n');
    fprintf (fid, '\\begin{axis}[\n');
    fprintf (fid, '    width=10cm, height=7cm,\n');
    fprintf (fid, '    xlabel={$u_r$ at outer boundary [mm]},\n');
    fprintf (fid, '    ylabel={Internal pressure $P$ [MPa]},\n');
    fprintf (fid, '    title={Eighth-sphere: load--displacement curve},\n');
    fprintf (fid, '    legend pos=north west, grid=both]\n');
    fprintf (fid, '\\addplot[black, thick, no marks] table[x=u, y=P] {sphere_load_disp_hill.dat};\n');
    fprintf (fid, '\\addlegendentry{Hill (analytical)}\n');
    for k = 1:numel(proj_methods)
        fprintf (fid, '\\addplot[%s] table[x=u, y=P] {sphere_load_disp_%s.dat};\n', ...
                 method_pgf_styles{k}, proj_methods{k});
        fprintf (fid, '\\addlegendentry{%s}\n', latex_escape(proj_methods{k}));
    end
    fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');

    % --- Radial-stress plots --------------------------------------------------
    for jload = 1:numel(load_steps)
        step = load_steps(jload);
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{axis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={Radial coordinate $r$ [mm]},\n');
        fprintf (fid, '    ylabel={$\\sigma_r$ [MPa]},\n');
        fprintf (fid, '    title={Radial stress, load step %d},\n', step);
        fprintf (fid, '    legend pos=south east, grid=both]\n');
        fprintf (fid, '\\addplot[black, thick, no marks] table[x=r, y=hill] {sphere_sigma_r_step%d.dat};\n', step);
        fprintf (fid, '\\addlegendentry{Hill}\n');
        for k = 1:numel(proj_methods)
            fprintf (fid, '\\addplot[%s, mark repeat=10] table[x=r, y=%s] {sphere_sigma_r_step%d.dat};\n', ...
                     method_pgf_styles{k}, proj_methods{k}, step);
            fprintf (fid, '\\addlegendentry{%s}\n', latex_escape(proj_methods{k}));
        end
        fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    % --- Tangential-stress plots ---------------------------------------------
    for jload = 1:numel(load_steps)
        step = load_steps(jload);
        fprintf (fid, '\\begin{tikzpicture}\n');
        fprintf (fid, '\\begin{axis}[\n');
        fprintf (fid, '    width=10cm, height=7cm,\n');
        fprintf (fid, '    xlabel={Radial coordinate $r$ [mm]},\n');
        fprintf (fid, '    ylabel={$\\sigma_t$ [MPa]},\n');
        fprintf (fid, '    title={Tangential stress, load step %d},\n', step);
        fprintf (fid, '    legend pos=south east, grid=both]\n');
        fprintf (fid, '\\addplot[black, thick, no marks] table[x=r, y=hill] {sphere_sigma_t_step%d.dat};\n', step);
        fprintf (fid, '\\addlegendentry{Hill}\n');
        for k = 1:numel(proj_methods)
            fprintf (fid, '\\addplot[%s, mark repeat=10] table[x=r, y=%s] {sphere_sigma_t_step%d.dat};\n', ...
                     method_pgf_styles{k}, proj_methods{k}, step);
            fprintf (fid, '\\addlegendentry{%s}\n', latex_escape(proj_methods{k}));
        end
        fprintf (fid, '\\end{axis}\n\\end{tikzpicture}\n\n');
    end

    fprintf (fid, '\\end{document}\n');
end

function s = latex_escape(s)
    s = strrep(s, '_', '\_');
end

function [errl2_rad, errl2_tan] = hspace_l2_error_stress (hspace, hmsh, sigma, sigma_ex)

errl2_rad = 0;
errl2_tan = 0;
last_dof = cumsum (hspace.ndof_per_level);

for ilev = 1:hmsh.nlevels
  if (hmsh.nel_per_level(ilev) > 0)
      msh_level = hmsh.msh_lev{ilev};
      sp_level = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
      
      [errl2_rad_lev, errl2_tan_lev] = ...
      sp_l2_error_stress (sp_level, msh_level, hspace.Csub{ilev}*sigma(1:last_dof(ilev),:), sigma_ex);

      errl2_rad = errl2_rad + errl2_rad_lev.^2;
      errl2_tan = errl2_tan + errl2_tan_lev.^2;
  end
end
errl2_rad  = sqrt (errl2_rad);
errl2_tan  = sqrt (errl2_tan);

end

function [errl2_rad, errl2_tan] = sp_l2_error_stress (sp, msh, sigma, sigma_ex)
  
  voigt   = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];
 
  for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
  end
  [sigma_rad_ex, sigma_tan_ex] = feval (sigma_ex, x{:});
  sigma_rad_ex  = reshape (sigma_rad_ex, sp.ncomp, msh.nqn, msh.nel);
  sigma_tan_ex  = reshape (sigma_tan_ex, sp.ncomp, msh.nqn, msh.nel);

  w = msh.quad_weights .* msh.jacdet;

  sigma_matrix = zeros(3,3,msh.nqn*msh.nel);

  for icomp = 1:6
      
      eu = sp_eval_msh (sigma(:,icomp), sp, msh);
      eu = reshape (eu, sp.ncomp, msh.nqn * msh.nel);

      sigma_matrix(voigt(icomp,1), voigt(icomp,2), :) = eu;
      if icomp > 3
          sigma_matrix(voigt(icomp,2), voigt(icomp,1), :) = eu;
      end
  end
  % sigma_matrix = permute(sigma_matrix,[1 3 2]);

  for iGP = 1:msh.nqn*msh.nel
    
    coordX = x{1}(iGP);
    coordY = x{2}(iGP);
    coordZ = x{3}(iGP);

    rad_dir = [coordX; coordY; coordZ] ./ norm([coordX; coordY; coordZ]);
    tan_dir = [coordY; -coordX; 0] ./ norm([coordX; coordY; 0]);
  
    sigma_rad(iGP) = rad_dir' * sigma_matrix(:,:,iGP) * rad_dir;
    sigma_tan(iGP) = tan_dir' * sigma_matrix(:,:,iGP)  * tan_dir;
  end

  sigma_rad = reshape(sigma_rad, [1, msh.nqn, msh.nel]);
  sigma_tan = reshape(sigma_tan, [1, msh.nqn, msh.nel]);
  
  errl2_rad_elem = sum (reshape (sum ((sigma_rad - sigma_rad_ex).^2, 1), [msh.nqn, msh.nel]) .* w);
  errl2_tan_elem = sum (reshape (sum ((sigma_tan - sigma_tan_ex).^2, 1), [msh.nqn, msh.nel]) .* w);

  errl2_rad  = sqrt (sum (errl2_rad_elem));
  errl2_tan  = sqrt (sum (errl2_tan_elem));
  
end

function [sigma_r, sigma_t] = eval_Hill_sphere(s_y, x, y, z, a, b, P)

P_0 = 2*s_y/3*(1-(a^3/b^3));
radius = sqrt(x.^2 + y.^2 + z.^2);
sigma_r =  zeros(size(radius));
sigma_t =  zeros(size(radius));

if P < P_0
    for j = 1:size(radius, 2)
        for i = 1:size(radius, 1)
            sigma_r(i,j) = -P*a^3 / (b^3 - a^3) * (b^3 / radius(i,j)^3 - 1);
            sigma_t(i,j) =  P*a^3 / (b^3 - a^3) * (0.5*b^3 / radius(i,j)^3 + 1);
        end
    end

else
    fun_front = @(x) -P +2*s_y* log(x/a) + 2/3*s_y* (1- x.^3/b^3);
    c = fsolve(fun_front, 100);

    for j =1:size(radius, 2)
        for i = 1:size(radius, 1)
            if radius(i,j)<=c
                sigma_r(i,j) = -2*s_y*( log(c/radius(i,j)) + 1/3* (1-c^3/b^3));
                sigma_t(i,j) = 2*s_y*( 0.5- log(c/radius(i,j)) - 1/3* (1-c^3/b^3));
            else
                sigma_r(i,j) = -2*s_y*c^3 / ( 3 * b^3) * (b^3 / radius(i,j)^3 -1);
                sigma_t(i,j) = 2*s_y*c^3 / ( 3 * b^3) * (0.5* b^3 / radius(i,j)^3 +1);
            end
        end
    end

end
end