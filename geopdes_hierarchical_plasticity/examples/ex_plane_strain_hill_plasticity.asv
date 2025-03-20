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

% Physical parameters
E  =  210000;                                  % MPa
nu = 0.3;                                      % -
problem_data.yield_stress = @(x, y) 240 * ones (size (x));  % MPa
problem_data.kappa_lame = @(x, y) E/(3*(1-2*nu)) * ones (size (x));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
P = 192.09*.98;                                        % Limit internal pressure [MPa]
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) test_plane_strain_ring_plasticity_g_nmnn (x, y, P, nu, ind);
problem_data.h = @(x, y, ind) test_plane_strain_ring_plasticity_uex (x, y, E, nu, P);
problem_data.p = @(x, y, ind) P * ones (size (x));

% Plot in Matlab Hill solution
figure(1)
[u_ex,P_ex] = Hill_solution (E, nu, problem_data.yield_stress(1), 100, 200);
plot (u_ex,P_ex,'-k');

xlabel('u');
ylabel('P');
hold on
grid on
drawnow

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
p = 2;
method_data.degree      = [p p];                          % Degree of the splines
method_data.regularity  = method_data.degree - 1;         % Regularity of the splines
method_data.nsub_coarse = [2 2];                          % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];                          % Number of subdivisions for each refinement
method_data.nquad       = [p+1 p+1];                      % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';                     % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;                              % 0: False, 1: True
method_data.nload      = 10;                              % Number of load steps
method_data.newton_tol = 1e-8;                            % Newton tolerance
method_data.newton_iter_max = 100;                        % Newton max number of iterations

adaptivity_data.flag = 'elements';
% adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'GR';
adaptivity_data.max_level = 3;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 8;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5;

% 3) CALL TO THE SOLVER
[geometry, cell_hmsh, cell_hspace, cell_u] = adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);



%% 4) POST-PROCESSING.
u_r=zeros(1,method_data.nload+1);
u_ex=zeros(1,method_data.nload+1);
P_i=zeros(1,method_data.nload+1);


for i=1:method_data.nload+1
    % Exact solution (optional)
    P_i(i) = P/method_data.nload*(i-1);
    if i > 1
        [eu, F] = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {100,50});
        u_r(i) = norm(eu);
    else
        u_r(i) = 0;
    end
    
end

% 4.1) Plot in Matlab numerical results
plot (u_r,P_i,'-x');
drawnow

% 4.2) Export to Paraview
output_file = strcat('plane_strain_ring_hier_',num2str(method_data.nload));
vtk_pts = {linspace(0, 100, 21), linspace(0, 100, 21)};

fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (cell_u{method_data.nload}, cell_hspace{method_data.nload}, geometry, vtk_pts, output_file, {'displacement'}, {'value'})