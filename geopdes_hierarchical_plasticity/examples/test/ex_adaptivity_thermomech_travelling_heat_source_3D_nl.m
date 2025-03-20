% PHYSICAL DATA OF THE PROBLEM
close all
clear problem_data_th problem_data_mec
clc
addpath(genpath('/home/blackmamba/workspace/GeoPDEs'));

% Physical domain, defined as NURBS map given in a text file
problem_data_th.geo_name = 'geo_square_3D.txt';
% Generate Output folder
problem_output.folder = '/home/blackmamba/workspace/GeoPDEs/output_thermomech_coarse_time_3D_not_admissible.02param_withmech';
mkdir(problem_output.folder);

% Set non-linear flag
% if true the non-linear solver is used
problem_data_th.flag_nl = true;
% Set lumped matrix
problem_data_th.lumped = false;

% Non-linear analysis
problem_data_th.Newton_tol = 1.0e-05;
problem_data_th.num_Newton_iter = 20;
problem_data_th.non_linear_convergence_flag = 1;

% Type of boundary conditions for each side of the domain
problem_data_th.nmnn_sides   = [];
problem_data_th.drchlt_sides = [];
problem_data_th.convection_sides = [];
problem_data_th.radiation_sides = [];

% Physical parameters
problem_data_th.c_diff  = @conductivity; %@(x,y,z) ones(size(x))*29;         %
problem_data_th.grad_c_diff =   @conductivity_der_3D; % @(x,y,z) cat (1, ...
           % reshape (zeros(size(x)), [1, size(x)]), ...
           % reshape (zeros(size(x)), [1, size(x)]), ...
           % reshape (zeros(size(x)), [1, size(x)]));
                                                            %
problem_data_th.c_cap = @capacity; %  @(x,y,z) ones(size(x))*7820*600;     % 
problem_data_th.initial_temperature = 80.0;                      %[°C]

problem_data_th.coeff_th_exp = 2.0E-5; %mm/(mm°C)

% Time discretization
n_time_steps = 100;%200;%
time_end = 1e-2;%2e-2; %[s]
problem_data_th.time_discretization = linspace(0.0, time_end, n_time_steps+1);

% Heat Source path
P = 180; %[W]
laser_radius = 50e-03; %[mm]
laser_penetration_depth = 15e-03; %[mm]
absorbtion_coeff = 0.35;
problem_data_th.x_begin = 2.0; %[mm]
problem_data_th.x_end = 3.0; %[mm]
num_hatch=5; %10;
x_path_1 = linspace(problem_data_th.x_begin, problem_data_th.x_end, n_time_steps/ num_hatch);
x_path_2 = linspace(problem_data_th.x_end, problem_data_th.x_begin, n_time_steps/ num_hatch);
% x_path = [x_path_1, x_path_2, x_path_1, x_path_2, x_path_1, x_path_2, x_path_1, x_path_2, x_path_1, x_path_2];
x_path = [x_path_1, x_path_2, x_path_1, x_path_2, x_path_1];

problem_data_th.y_begin = 2.4;
problem_data_th.y_end = 2.6;
hatch_space = laser_radius;
y_path_1 = linspace(problem_data_th.y_begin, problem_data_th.y_end, n_time_steps / num_hatch);
y_path_2 = linspace(problem_data_th.y_begin + hatch_space, problem_data_th.y_end +hatch_space, n_time_steps / num_hatch);
y_path_3 = linspace(problem_data_th.y_begin + hatch_space * 2, problem_data_th.y_end +hatch_space * 2, n_time_steps / num_hatch);
y_path_4 = linspace(problem_data_th.y_begin + hatch_space * 3, problem_data_th.y_end +hatch_space * 3, n_time_steps / num_hatch);
y_path_5 = linspace(problem_data_th.y_begin + hatch_space * 4, problem_data_th.y_end +hatch_space * 4, n_time_steps / num_hatch);
% y_path_6 = linspace(problem_data.y_begin + hatch_space * 5, problem_data.y_end +hatch_space * 5, n_time_steps / num_hatch);
% y_path_7 = linspace(problem_data.y_begin + hatch_space * 6, problem_data.y_end +hatch_space * 6, n_time_steps / num_hatch);
% y_path_8 = linspace(problem_data.y_begin + hatch_space * 7, problem_data.y_end +hatch_space * 7, n_time_steps / num_hatch);
% y_path_9 = linspace(problem_data.y_begin + hatch_space * 8, problem_data.y_end +hatch_space * 8, n_time_steps / num_hatch);
% y_path_10 = linspace(problem_data.y_begin + hatch_space * 9, problem_data.y_end +hatch_space * 9, n_time_steps / num_hatch);

% y_path = [y_path_1, y_path_2, y_path_3, y_path_4, y_path_5, y_path_6, y_path_7, y_path_8, y_path_9, y_path_10];
y_path = [y_path_1, y_path_2, y_path_3, y_path_4, y_path_5];

problem_data_th.z_begin = 1;
problem_data_th.z_end = 1;
z_path = linspace(problem_data_th.z_begin, problem_data_th.z_end, n_time_steps);
problem_data_th.path = [x_path', y_path', z_path'];

% Radiation and Convection Boundaries
problem_data_th.rad =  0.0;
problem_data_th.u_r = 3000.0;
problem_data_th.rad_fun = @(x,y,z,u,ind) problem_data_th.rad*(problem_data_th.u_r*ones(size(u))- u);
problem_data_th.rad_tilda =  0.0;
problem_data_th.rad_tilda_fun =  @(x,y,z,ind) ones(size(x)) * problem_data_th.rad_tilda;
problem_data_th.conv =  0.0;
problem_data_th.conv_fun = @(x,y,z,ind) ones(size(x)) * problem_data_th.conv;
problem_data_th.u_e = 20.0;
problem_data_th.conv_fun_rhs = @(x,y,z,u,ind) problem_data_th.conv*(problem_data_th.u_e*ones(size(u))- u);
problem_data_th.alpha = 1;  % alpha = 0 explicit Euler; alpha=1/2 Crank-Nicholson; aplha=1 backward Euler

% Source and boundary terms
problem_data_th.f = @(x,y,z,path_x,path_y,path_z) absorbtion_coeff * 6.0*sqrt(3.0)*P/(pi * sqrt(pi) * laser_radius * laser_radius * laser_penetration_depth) * ...
    exp(- 3*(x-path_x).^2/laser_radius^2 - 3*(y-path_y).^2/laser_radius^2 - 3*(z-path_z).^2/(laser_penetration_depth)^2);         % Body Load
problem_data_th.h = @(x, y, z, ind) ones(size(x))*problem_data_th.initial_temperature;                                                  % Dirichlet Boundaries

%% =============== MECHANICAL PROBLEM =====================================
% Physical domain, defined as NURBS map given in a text file
problem_data_mec.geo_name = 'geo_square_3D.txt';
problem_data_mec.flag_nl = false;
problem_data_mec.non_linear_convergence_flag = 0;

% Type of boundary conditions for each side of the domain
problem_data_mec.nmnn_sides   = [];
problem_data_mec.press_sides  = [];
problem_data_mec.drchlt_sides = 5;
problem_data_mec.symm_sides   = [];

% Physical parameters
E  =  195E+3; nu = 0.3; 
problem_data_mec.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
problem_data_mec.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Time discretization
n_time_steps = 100;
time_end = 1e-2; %[s]
problem_data_mec.time_discretization = linspace(0.0, time_end, n_time_steps+1);


% Source and boundary terms
fx = @(x, y, z) zeros(size(x));
fy = @(x, y, z) zeros(size(y));
fz = @(x, y, z) zeros(size(z));
problem_data_mec.f = @(x, y, z) cat(1, ...
                    reshape (fx (x,y,z), [1, size(x)]), ...
                    reshape (fy (x,y,z), [1, size(x)]), ...
                    reshape (fz (x,y,z), [1, size(x)]));
problem_data_mec.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data_th
method_data_th.degree      = [2 2 2];        % Degree of the splines
method_data_th.regularity  = [1 1 1];        % Regularity of the splines
method_data_th.nsub_coarse = [1 1 1];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data_th.nsub_refine = [2 2 2];        % Number of subdivisions for each refinement
method_data_th.nquad       = [3 3 3];        % Points for the Gaussian quadrature rule
method_data_th.space_type  = 'standard';     % 'simplified' (only children functions) or 'standard' (full basis)
method_data_th.truncated   = 1;              % 0: False, 1: True
clear method_data_mec
method_data_mec.degree      = [2 2 2];        % Degree of the splines
method_data_mec.regularity  = [1 1 1];        % Regularity of the splines
method_data_mec.nsub_coarse = [1 1 1];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data_mec.nsub_refine = [2 2 2];        % Number of subdivisions for each refinement
method_data_mec.nquad       = [3 3 3];        % Points for the Gaussian quadrature rule
method_data_mec.space_type  = 'standard';     % 'simplified' (only children functions) or 'standard' (full basis)
method_data_mec.truncated   = 1;              % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data_th
adaptivity_data_th.flag = 'elements';
adaptivity_data_th.doCoarsening = true;
adaptivity_data_th.C0_est = 1.0;
adaptivity_data_th.mark_neighbours = true;
adaptivity_data_th.mark_param = 0.02;% 0.5; %
adaptivity_data_th.mark_param_coarsening = 0.1;
adaptivity_data_th.crp = 2.0;                     %coarsening relaxation parameter
adaptivity_data_th.radius = [laser_radius, laser_radius, laser_penetration_depth / 2]*2.0;
adaptivity_data_th.strategy = 'adapt'; % adapt = error-based adaptivity or 'geom' = geometric based adaptivity
adaptivity_data_th.adm_strategy = 'admissible'; % 'admissible' or 'balancing'
adaptivity_data_th.adm = 9; %10; % 1 + method_data.truncated;%  
adaptivity_data_th.coarsening_flag = 'any'; %'any', 'all'
adaptivity_data_th.coarse_flag = 'L2_global'; % 'bezier', 'MS_all', 'MS_old', 'L2_global'
adaptivity_data_th.mark_strategy = 'MS';
adaptivity_data_th.max_level = 9;
adaptivity_data_th.max_ndof = 1000000;
adaptivity_data_th.num_max_iter = 12;
adaptivity_data_th.max_nel = 100000;
adaptivity_data_th.tol = 1.0e-04;
adaptivity_data_th.timeToRefine = linspace(1,n_time_steps,n_time_steps);

clear adaptivity_data_mec
adaptivity_data_mec.flag = 'elements';
adaptivity_data_mec.doCoarsening = false;
adaptivity_data_mec.C0_est = 1.0;
adaptivity_data_mec.mark_neighbours = true;
adaptivity_data_mec.mark_param = 0.75;% 0.5; %
adaptivity_data_mec.mark_param_coarsening = 0.1;
adaptivity_data_mec.crp = 2.0;                     %coarsening relaxation parameter
adaptivity_data_mec.radius = [laser_radius, laser_radius, laser_penetration_depth / 2]*2.0;
adaptivity_data_mec.strategy = 'adapt'; % adapt = error-based adaptivity or 'geom' = geometric based adaptivity
adaptivity_data_mec.adm_strategy = 'admissible'; % 'admissible' or 'balancing'
adaptivity_data_mec.adm = 2; %10; % 1 + method_data.truncated;%  
adaptivity_data_mec.coarsening_flag = 'any'; %'any', 'all'
adaptivity_data_mec.coarse_flag = 'L2_global'; % 'bezier', 'MS_all', 'MS_old', 'L2_global'
adaptivity_data_mec.mark_strategy = 'GR';
adaptivity_data_mec.max_level = 6;
adaptivity_data_mec.max_ndof = 10000000;
adaptivity_data_mec.num_max_iter = 5;
adaptivity_data_mec.max_nel = 1000000;
adaptivity_data_mec.tol = 1.0e-04;
adaptivity_data_mec.timeToRefine = 1;
adaptivity_data_mec.timeToSolveMec = [1,11,21,31,41,51,61,71,81,91,100];

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.adaptivity = true;
plot_data.print_info = true;
plot_data.plot_matlab = false;
plot_data.time_steps_to_post_process = [1,2]; %linspace(0,n_time_steps,21);%[1,2,3,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,...
  %  525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,999,1000];%linspace(1,n_time_steps,n_time_steps);%[1,10,25,50,75,100,150,200,250,300,350,400];%200; % [1,25,50,75,100,150,200,250,300,350,400];% [1,5,10,15,20];
plot_data.file_name = strcat(problem_output.folder, '/test_moving_heat_source_thermomech_3D_%d');
plot_data.file_name_mesh = strcat(problem_output.folder, '/test_moving_heat_source_thermomech_3D_mesh_%d/hmsh');
plot_data.folder_name_mesh = strcat(problem_output.folder, '/test_moving_heat_source_thermomech_3D_mesh_%d');
plot_data.file_name_dofs = strcat(problem_output.folder, '/test_moving_heat_source_thermomech_3D_var.csv');
plot_data.file_name_temp_plot = strcat(problem_output.folder, '/test_moving_heat_source_thermomech_3D_temperature_%d');
plot_data.file_name_var = strcat(problem_output.folder, '/test_moving_heat_source_thermal_3D_variables_%d.mat');
plot_data.file_name_varMec = strcat(problem_output.folder, '/test_moving_heat_source_mechanical_3D_variables_%d.mat');

plot_data.npoints_x = 51;        %number of points x-direction in post-processing
plot_data.npoints_y = 51;        %number of points x-direction in post-processing
plot_data.npoints_z = 11;        %number of points x-direction in post-processing

fid = fopen (plot_data.file_name, 'w');
tic
[geometry, hmsh, hspace, u, solution_data] = adaptivity_thermomech_transient(problem_data_th, problem_data_mec, method_data_th, method_data_mec, adaptivity_data_th, adaptivity_data_mec, plot_data);
toc
fprintf (fid, num2str(toc));
fclose(fid);
