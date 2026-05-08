% Approximating a function f:\RR^2\rightarrow \RR^m with the Quasi-Interpolant

% Here load hierarchical space, hierarchical mesh and data, if necessary
% data must be an Nx(2+m) matrix
% The first two columns of data are points in the parametric domain
% the remaining comlumns are the values of the approximated vector-valued function

% The parametric points (first two columns of data) 
% need to be re-scaled so that they fall into [0,1]x[0,1]
a=min(data(:,1));
b=max(data(:,1));
c=min(data(:,2));
d=max(data(:,2));
data(:,1)=(data(:,1)-a)/(b-a);
data(:,2)=(data(:,2)-c)/(d-c);

% Computing the coefficients of the quasi-interpolant in the space
% hspace based on the mesh hmsh, approximating the data.
% libqi backend (C/OpenMP). See ./libqi/README.md for build instructions
% (run libqi/build_mex.m once to compile qi_mex for your platform).
addpath(fullfile(fileparts(mfilename('fullpath')), 'libqi', 'matlab'));
QI_coeff = get_QI_coeffs_c(hspace, hmsh, ...
                           struct('x', data(:,1), ...
                                  'y', data(:,2), ...
                                  'f', data(:,3:end)));

%evaluation on scattered data...
QI_val_scattered=sp_eval_alt (QI_coeff, hspace, data(:,1:2));
%...or evaluation on mxn grid 
m=100; n=100; 
for k=1:size(QI_coeff,2) %cycle over the number of components of the approximated function
    QI_val_grid(k,:,:)=sp_eval (QI_coeff(:,k), hspace, geometry,[m n]); 
end
