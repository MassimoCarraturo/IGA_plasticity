% Clear QI_graded results from checkpoint so they re-run with new L2 projection
cp = fullfile(fileparts(mfilename('fullpath')), 'results', 'qi_study_v2', 'qi_study_checkpoint.mat');
tmp = load(cp, 'R');
R = tmp.R;
tags = {'QI_graded_lev1', 'QI_graded_lev2', 'QI_graded_lev3', 'QI_graded_lev4'};
for i = 1:numel(tags)
    if isfield(R, tags{i})
        R.(tags{i}).ndof = 0;
        fprintf('  Cleared %s\n', tags{i});
    end
end
save(cp, 'R');
fprintf('Checkpoint updated.\n');
