% EX_PLANE_STRAIN_ADAPTIVITY: solve the plane-strain problem on one quarter of a cylinder.

% PHYSICAL DATA OF THE PROBLEM
clear all
clc
close all
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_ring_SouzaNeto.txt';

% Type of boundary conditions
problem_data.nmnn_sides   = [1];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [];
problem_data.symm_sides   = [3 4];
problem_data.slider_sides = [];

% Physical parameters
E  =  210000;                                  % MPa
nu = 0.3;                                      % -
sigma_y = 240;
mu_lame =  E/(2*(1+nu));
kappa_lame = E/(3*(1-2*nu));
problem_data.yield_stress = @(x, y) sigma_y * ones (size (x));  % MPa
problem_data.kappa_lame = @(x, y) kappa_lame * ones (size (x));
problem_data.mu_lame = @(x, y) mu_lame * ones (size (x));

% Source and boundary terms
P = 192.09*.98;                                        % Limit internal pressure [MPa]
nload = 10;
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) test_plane_strain_ring_plasticity_g_nmnn (x, y, P, nu, ind);
problem_data.h = @(x, y, ind) test_plane_strain_ring_plasticity_uex (x, y, E, nu, P);
problem_data.p = @(x, y, ind) P * ones (size (x));

% Plot in Matlab Hill solution

[u_ex,P_ex, sigma_r_ex, sigma_t_ex, radius] = Hill_solution (E, nu, problem_data.yield_stress(1), 100, 200,P, nload);
load_step_eval_stress = [1, 10]; %[nload, int64(nload/2)]*int64(100/nload);
figure(1)
plot (u_ex,P_ex,'-k');
xlabel('u');
ylabel('P');
hold on
grid on
drawnow

figure(2)
xlabel('radial coord');
ylabel('sigma radial');
hold on
for jload =1:length(load_step_eval_stress)
    disp(jload)
    plot (radius,sigma_r_ex(load_step_eval_stress(jload),:),'-k');    
end
grid on
drawnow

figure(3)
xlabel('radial coord');
ylabel('sigma tangential');
hold on
for jload =1:length(load_step_eval_stress)
    plot (radius,sigma_t_ex(load_step_eval_stress(jload),:),'-k');    
end
grid on
drawnow




% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
p = 3;
method_data.degree      = [p p];                          % Degree of the splines
method_data.regularity  = method_data.degree - 1;         % Regularity of the splines
method_data.nsub_coarse = [4 4];                          % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];                          % Number of subdivisions for each refinement
method_data.nquad       = [p+1 p+1];                      % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';                     % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;                              % 0: False, 1: True
method_data.nload      = nload;                              % Number of load steps
method_data.newton_tol = 1e-8;                            % Newton tolerance
method_data.newton_tol_abs = 1e-10;                            % Newton tolerance
method_data.newton_iter_max = 100;                        % Newton max number of iterations
method_data.type_projection = 'QI';                         % 'QI' / 'L2'

adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .8;
adaptivity_data.mark_param_coarsening = .05;
adaptivity_data.mark_strategy = 'MS'; % GR/MS/GERS
adaptivity_data.max_level = 4;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5 *3.14 *30000;

adaptivity_data.adm_strategy = 'admissible';
adaptivity_data.coarsening_flag = 'any'; %'any', 'all'
adaptivity_data.adm = p ; %1 + method_data.truncated;
% adaptivity_data.coarse_flag = 'bezier'; % 'bezier', 'MS_all', 'MS_old', 'L2_global'
% 3) CALL TO THE SOLVER
[geometry, cell_hmsh, cell_hspace,  cell_hspace_scalar,  cell_u, cell_eps_pl, cell_sigma, solution_data] = adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);



%% 4) POST-PROCESSING.


% displacements
u_r=zeros(1,method_data.nload+1);
P_i=zeros(1,method_data.nload+1);
for i=1:method_data.nload+1
    P_i(i) = P/method_data.nload*(i-1);
    if i > 1
        [eu, F] = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1,.5});
        u_r(i) = norm(eu);
    else
        u_r(i) = 0;
    end    
end


% stress
n_points_radial = 100;
pt_eval ={linspace(0,1,n_points_radial),0.5};
sigma_r = zeros(length(load_step_eval_stress), n_points_radial);
sigma_t = zeros(length(load_step_eval_stress), n_points_radial);

for jload =1:length(load_step_eval_stress)
    eps_tot = zeros(6, n_points_radial);
    eps_pl = zeros(6, n_points_radial);
    sigma = zeros(6, n_points_radial);
    

    % compute epsilon total at evaluation points
    [eu, F] = sp_eval (cell_u{load_step_eval_stress(jload)}, cell_hspace{load_step_eval_stress(jload)}, geometry, pt_eval  , 'gradient');
    eu = reshape(eu, [2,2,n_points_radial]);
    eps_tot_plane = .5*(eu +permute(eu,[2,1,3]));
    eps_tot(1,:) = eps_tot_plane(1,1,:);
    eps_tot(2,:) = eps_tot_plane(2,2,:);
    eps_tot(4,:) = eps_tot_plane(1,2,:);

    % compute epsilon plastice at evaluation points
    for comp =1:6
        [eu, F] = sp_eval (cell_eps_pl{load_step_eval_stress(jload)}(:,comp), cell_hspace_scalar{load_step_eval_stress(jload)}, geometry, pt_eval );
        eps_pl(comp,:)=eu;
    end
    
    % compute sigma
    for ipt = 1:n_points_radial
        [SIGMA, ~, ~] = perfect_plasticity_model(eps_tot(:,ipt), eps_pl(:,ipt), mu_lame, kappa_lame, sigma_y);
        sigma(:,ipt) = SIGMA;
        sigma_r(jload,ipt) = .5*(sigma(1,ipt)+sigma(2,ipt) + 2*sigma(4,ipt));
        sigma_t(jload,ipt) = .5*(sigma(1,ipt)+sigma(2,ipt) - 2*sigma(4,ipt));
    end

    %

end


% 4.1) Plot in Matlab numerical results
figure(1)
plot (u_r,P_i,'-x');
drawnow

figure(2)
for jload =1:length(load_step_eval_stress)
    plot (linspace(100,200,n_points_radial),sigma_r(jload,:),'-x');
end
drawnow

figure(3)
for jload =1:length(load_step_eval_stress)
    plot (linspace(100,200,n_points_radial),sigma_t(jload,:),'-x');
end
drawnow

% 4.2) Export to Paraview
output_file = strcat('plane_strain_ring_hier_',method_data.type_projection,'_',num2str(method_data.nload));
vtk_pts = {linspace(0, 1, 41), linspace(0, 1, 41)};

fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (cell_u{method_data.nload}, cell_hspace{method_data.nload}, geometry, vtk_pts, output_file, {'displacement'}, {'value'})



hspace_scalar = cell_hspace_scalar{method_data.nload};
for comp =1:6
    output_file = strcat('plane_strain_ring_hier_',method_data.type_projection,'_',num2str(method_data.nload),'_eps_pl_',num2str(comp) );    
    name_field = strcat('plastic_strain_',num2str(comp));
    sp_to_vtk ( cell_eps_pl{method_data.nload}(:,comp), hspace_scalar, geometry, vtk_pts, output_file, {name_field}, {'value'})

    output_file = strcat('plane_strain_ring_hier_',method_data.type_projection,'_',num2str(method_data.nload),'_sigma_',num2str(comp) );    
    name_field = strcat('sigma_',num2str(comp));
    sp_to_vtk ( cell_sigma{method_data.nload}(:,comp), hspace_scalar, geometry, vtk_pts, output_file, {name_field}, {'value'})
end




