function results_dir = fm_setup (here, study)
% FM_SETUP  Paths, warnings and results directory of the front-metrics study drivers.
% The FM_DEGREE, FM_NU and FM_TRANSFER overrides (see fm_env) each add a suffix to the directory name
  project_root = fullfile (here, '..', '..', '..');
  addpath(genpath(fullfile(project_root,'nurbs-1.4.3','nurbs-1.4.3','inst')));
  addpath(genpath(fullfile(project_root,'geopdes-3.2.2','geopdes','inst')));
  addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','hierarchical_classes')));
  addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','adaptivity_iga')));
  addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','initialize')));
  addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','examples','plasticity')));
  addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','quasi_interpolant_hierarchical')));
  addpath(fullfile(project_root,'geopdes_hierarchical_plasticity','quasi_interpolant_hierarchical','libqi','matlab'));
  warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');
  results_dir = fullfile (here, 'results', [study fm_env('DEGREE', '', @(s) ['_p' s]) ...
                fm_env('NU', '', @(s) ['_nu' strrep(s, '.', '')]) fm_env('TRANSFER', '', @(s) ['_' s])]);
  if (~exist (results_dir, 'dir')); mkdir (results_dir); end
end
