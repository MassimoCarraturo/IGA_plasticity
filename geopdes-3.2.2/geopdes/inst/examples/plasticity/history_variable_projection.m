function eps_pl_control_var = history_variable_projection(space, msh, eps_pl, type_proj)

n_hist_var = size(eps_pl, 3);
eps_pl_control_var= zeros(space.ndof, n_hist_var);

if  strcmpi(type_proj, 'L2')
    % global L2 projection
    M = op_u_v_tp(space, space, msh);

    for ivar = 1:n_hist_var
        % compute RHS
        for idim = 1:msh.ndim
            size1 = size (space.sp_univ(idim).connectivity);
            if (size1(2) ~= msh.nel_dir(idim))
                error ('The discrete space is not associated to the mesh')
            end
        end
        rhs = zeros (space.ndof, 1);
        for iel = 1:msh.nel_dir(1)
            msh_col = msh_evaluate_col (msh, iel);
            sp_col  = sp_evaluate_col (space, msh_col);        
            rhs = rhs + op_f_v (sp_col, msh_col,  eps_pl(msh_col.elem_list,:,ivar)'  );
        end
        
        % solve L2
        eps_pl_control_var(:,ivar) = M\rhs;
    end

elseif  strcmpi(type_proj, 'QI')

    
    
    n_quad_nodes = size(eps_pl,1) *  size(eps_pl,2);
    data = zeros(n_quad_nodes, msh.ndim);
    f = zeros(n_quad_nodes, n_hist_var);

    
  
    % msh_val = msh_precompute(msh);
    quad_nodes = my_msh_evaluate_qn (msh, 1:msh.nel);
    data(:,1:msh.ndim ) = reshape(quad_nodes, msh.ndim, n_quad_nodes)';
    f(:, :) = reshape(permute(eps_pl, [3,2,1]),  n_hist_var , n_quad_nodes)';


    hmsh     = hierarchical_mesh (msh);
    hspace   = hierarchical_space (hmsh, space,  'standard', 1, space.degree-1);
    
    QI_coeff = getcoeff_localLS_Bspl(hspace,hmsh,data, f, 0.);
    eps_pl_control_var(:,:) = QI_coeff;


end



end