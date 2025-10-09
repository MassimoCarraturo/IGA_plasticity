% EX_PLANE_STRAIN_RING: Benchmark taken from "Computationa Methods for
% Plasticity: Theory and Applications" EA de Souza Neto, D. Reric, DRJ
% Owen, Eds. Wiley, pp.244-247

% 1) PHYSICAL DATA OF THE PROBLEM
clear all
clc
close all
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
problem_data.yield_stress = @(x, y, z) 240 * ones (size (x));  % MPa
problem_data.kappa_lame = @(x, y, z) E/(3*(1-2*nu)) * ones (size (x));
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
nload = 5;
P = 332*.99;                                        % Limit internal pressure [MPa]
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2), size (x, 3)   );
problem_data.p = @(x, y, z) P*ones (size (x));

% Plot in Matlab Hill solution
[u_ex,P_ex, sigma_r_ex, sigma_t_ex, radius] = sphere_solution (E, nu, problem_data.yield_stress(1), 100, 200,P,nload);
load_step_eval_stress = [2, 5];%[5,10]; %[nload, int64(nload/2)]*int64(100/nload);

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
for p =2:2 % loop for p-refinemet study
    clear method_data
    method_data.degree     = [p p p];                      % Degree of the basis functions
    method_data.regularity = [p-1 p-1 p-1];     % Regularity of the basis functions
    method_data.nsub_coarse = [15,5,5];            % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
    method_data.nsub_refine = [2 2 2];                          % Number of subdivisions for each refinement
    method_data.nquad      = [p+1 p+1 p+1];                  % Points for the Gaussian quadrature rule
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
    adaptivity_data.mark_param = 0.9;
    adaptivity_data.mark_param_coarsening = .05;
    adaptivity_data.mark_strategy = 'MS'; % GR/MS/GERS
    adaptivity_data.max_level = 1;
    adaptivity_data.max_ndof = 15000;
    adaptivity_data.num_max_iter = 1;
    adaptivity_data.max_nel = 5000;
    adaptivity_data.tol = 1e-5 *3.14 *30000;
    
    adaptivity_data.adm_strategy = 'admissible';
    adaptivity_data.coarsening_flag = 'any'; %'any', 'all'
    adaptivity_data.adm = p ;

    % 3) CALL TO THE SOLVER
    [geometry, cell_hmsh, cell_hspace,  cell_hspace_scalar,  cell_u, cell_eps_pl, cell_sigma, solution_data] = ...
    adaptivity_J2_plasticity (problem_data, method_data, adaptivity_data);

   %% 4) POST-PROCESSING.
   % displacements
    u_r=zeros(1,method_data.nload+1);
    P_i=zeros(1,method_data.nload+1);
    for i=1:method_data.nload+1
        P_i(i) = P/method_data.nload*(i-1);
        if i > 1
            [eu, F] = sp_eval (cell_u{i-1}, cell_hspace{i-1}, geometry, {1,.5,.5});
            u_r(i) = norm(eu);
        else
            u_r(i) = 0;
        end    
    end

    % 4.1) Plot in Matlab numerical results
    figure(1)
    plot (u_r,P_i,'-x');
    drawnow

    % compute radial and tangential stresses
    rad_dir = [.5;.5;sqrt(2)/2];
    tan_dir =  [.5;.5;-sqrt(2)/2];
    voigt = [1,1;2,2; 3,3; 1,2; 2,3;1,3];
    rad_pos = linspace(0,1,length(radius));

    sigma_rad = zeros(length(radius),length(load_step_eval_stress) );
    sigma_tan = zeros(length(radius),length(load_step_eval_stress) );

    for jload =1:length(load_step_eval_stress)
    for jpoint =1:length(rad_pos)
        position = {rad_pos(jpoint), 0.5, 0.75};        
        sigma_matrix = zeros(3,3);
        for icomp =1:6
            [eu, F] = sp_eval (cell_sigma{load_step_eval_stress(jload)}(:,icomp), cell_hspace_scalar{load_step_eval_stress(jload)}, geometry, position);
            sigma_matrix(voigt(icomp,1), voigt(icomp,2)) = eu; 
            if icomp > 3
                sigma_matrix(voigt(icomp,2), voigt(icomp,1)) = eu;   
            end
        end
        sigma_rad(jpoint, jload) = rad_dir' * sigma_matrix * rad_dir;
        sigma_tan(jpoint, jload) = tan_dir' * sigma_matrix * tan_dir;

    end
    end
    % plot stress
    figure(2)
    for jload =1:length(load_step_eval_stress)
        plot (radius,sigma_rad(:, jload),'--');    
    end

    figure(3)
    for jload =1:length(load_step_eval_stress)
        plot (radius,sigma_tan(:, jload),'--');    
    end


% 
% end % loop p refinement study
% legend ('Exact solution',  'p=2', 'p=3', 'p=4');
% 
% %error_l2 = sp_l2_error (space, msh, u(:,i), problem_data.uex)
end
%% 4.2) Export to Paraview
output_file = strcat('hier_solution_plastic_sphere_',num2str(i-1));
vtk_pts = {linspace(0, 1, 41), linspace(0, 1, 41), linspace(0, 1, 41)};

fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (cell_u{method_data.nload}, cell_hspace{method_data.nload}, geometry, vtk_pts, output_file, {'displacement'}, {'value'})



hspace_scalar = cell_hspace_scalar{method_data.nload};
for comp =1:6
    output_file = strcat('hier_solution_plastic_sphere_',method_data.type_projection,'_',num2str(method_data.nload),'_eps_pl_',num2str(comp) );    
    name_field = strcat('plastic_strain_',num2str(comp));
    sp_to_vtk ( cell_eps_pl{method_data.nload}(:,comp), hspace_scalar, geometry, vtk_pts, output_file, {name_field}, {'value'})

    output_file = strcat('hier_solution_plastic_sphere_',method_data.type_projection,'_',num2str(method_data.nload),'_sigma_',num2str(comp) );    
    name_field = strcat('sigma_',num2str(comp));
    sp_to_vtk ( cell_sigma{method_data.nload}(:,comp), hspace_scalar, geometry, vtk_pts, output_file, {name_field}, {'value'})
end
