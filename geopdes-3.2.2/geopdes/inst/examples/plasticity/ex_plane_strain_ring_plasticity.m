% EX_PLANE_STRAIN_RING: Benchmark taken from "Computationa Methods for
% Plasticity: Theory and Applications" EA de Souza Neto, D. Reric, DRJ
% Owen, Eds. Wiley, pp.244-247

% 1) PHYSICAL DATA OF THE PROBLEM
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
for p =1:4 % loop for p-refinemet study
    clear method_data
    method_data.degree     = [p p];                      % Degree of the basis functions
    method_data.regularity = method_data.degree - 1;     % Regularity of the basis functions
    method_data.nsub       = [10 10];                    % Number of subdivisions
    method_data.nquad      = [p+1 p+1];                  % Points for the Gaussian quadrature rule
    method_data.nload      = 10;                         % Number of load steps
    method_data.newton_tol = 1e-8;                       % Newton tolerance
    method_data.newton_iter_max = 100;                   % Newton max number of iterations

    % 3) CALL TO THE SOLVER
    [geometry, msh, space, u] = solve_J2_plasticity (problem_data, method_data);

    %% 4) POST-PROCESSING.
    u_r=zeros(1,method_data.nload+1);
    u_ex=zeros(1,method_data.nload+1);
    P_i=zeros(1,method_data.nload+1);
    u_plot = [zeros(size(u(:,1))) u];

    for i=1:method_data.nload+1
        % Exact solution (optional)
        P_i(i) = P/method_data.nload*(i-1);
        % problem_data.uex =  test_plane_strain_ring_plasticity_uex (x, y, E, nu, problem_data.yield_stress(1), P_i(i));
        [eu, F] = sp_eval (u_plot(:,i), space, geometry, {100,50});
        u_r(i) = norm(eu);
    end

    % 4.1) Plot in Matlab numerical results
    plot (u_r,P_i,'-x');
    drawnow

end % loop p refinement study
legend ('Exact solution', 'p=1', 'p=2', 'p=3', 'p=4');

%error_l2 = sp_l2_error (space, msh, u(:,i), problem_data.uex)

% 4.2) Export to Paraview
output_file = strcat('plane_strain_ring_Deg3_Reg2_Sub9_',num2str(i-1));
vtk_pts = {linspace(0, 100, 21), linspace(0, 100, 21)};

fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u_plot(:,i), space, geometry, vtk_pts, output_file, {'displacement'}, {'value'})