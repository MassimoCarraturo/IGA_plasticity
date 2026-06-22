% Wrapper to run the 2D cylinder QI-variant study with diary logging.
% NOTE: study_cylinder_qi_variants starts with `clear`, which wipes this
% workspace — so do not reference any local variable after the call.
out_dir = fullfile(fileparts(mfilename('fullpath')), 'results', 'cylinder_qi_study');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
diary(fullfile(out_dir, 'run_log.txt'));
fprintf('=== Cylinder QI study started ===\n');
try
    study_cylinder_qi_variants;
catch ME
    fprintf(2, 'ERROR: %s\n%s\n', ME.message, getReport(ME));
end
diary off;
