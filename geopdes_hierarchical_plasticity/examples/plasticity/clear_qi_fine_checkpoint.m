% Clear QI_fine entries from the cylinder checkpoint so they re-run
% (with the new finer-mesh regularization).  Other methods are kept.
cp = fullfile(fileparts(mfilename('fullpath')), 'results', 'cylinder_qi_study', 'cylinder_qi_checkpoint.mat');
tmp = load(cp, 'R');
R = tmp.R;
tags = {'QI_fine_lev1', 'QI_fine_lev2', 'QI_fine_lev3', 'QI_fine_lev4'};
for i = 1:numel(tags)
    if isfield(R, tags{i})
        R.(tags{i}).ndof = 0;
        fprintf('  Cleared %s\n', tags{i});
    end
end
save(cp, 'R');
fprintf('Checkpoint updated — QI_fine will re-run.\n');
