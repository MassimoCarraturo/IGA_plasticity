
function quad_nodes = my_msh_evaluate_qn (msh, elem_list)

elem_list = elem_list(:)';

msh_col.ndim = msh.ndim;
msh_col.rdim = msh.rdim;

msh_col.nel  = numel (elem_list);
msh_col.elem_list = elem_list;
msh_col.nel_dir = msh.nel_dir;

msh_col.nqn_dir = msh.nqn_dir;
msh_col.nqn  = msh.nqn;


indices = cell (msh.ndim, 1);
% Trailing 1 keeps ind2sub happy for 1D meshes on MATLAB R2024b+
[indices{:}] = ind2sub ([msh.nel_dir, 1], elem_list);
indices = cell2mat (indices);

qn_elems = arrayfun(@(ii) {msh.qn{ii}(:,indices(ii,:))}, 1:msh.ndim);
qqn = cell (1,msh_col.nel);
for iel = 1:numel(elem_list)
for idim = 1:msh.ndim
  qqn{iel}{idim} = qn_elems{idim}(:,iel)';
end
end

quad_nodes = zeros (msh.ndim, msh.nqn, numel(elem_list));
for iel = 1:numel(elem_list)
    xx = cell (msh.ndim, 1);
    [xx{:}] = ndgrid (qqn{iel}{:});
    for idim = 1:msh.ndim
      quad_nodes(idim,:,iel) = xx{idim}(:)';
    end
end

end