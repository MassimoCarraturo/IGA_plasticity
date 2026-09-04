% Wrapper to run the QI variant study with diary logging
cd(fileparts(mfilename('fullpath')));
logfile = fullfile(pwd, 'results', 'qi_study_v2', 'study_log.txt');
diary(logfile);
try
    study_sphere_qi_variants;
catch ME
    fprintf(2, 'ERROR: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf(2, '  in %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
end
diary off;
