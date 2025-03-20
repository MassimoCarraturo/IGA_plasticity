function eps_pl_control_var = history_variable_projection_hier(hspace, hmsh, eps_pl, type_proj)


n_levels = hmsh.nlevels;
n_hist_var = size(eps_pl{n_levels}, 3);
eps_pl_control_var= zeros(hspace.ndof, n_hist_var);

if  strcmpi(type_proj, 'L2')
    
    % not yet implemented

elseif  strcmpi(type_proj, 'QI')
    if hspace.ncomp == 3
        print('QI for 2D only!!!!!!!')
        return
    end
    

    n_quad_nodes = 0;
    for ilev = 1:n_levels
         if (hmsh.nel_per_level(ilev) > 0)
            n_quad_nodes = n_quad_nodes + hmsh.mesh_of_level(ilev).nqn * hmsh.nel_per_level(ilev);
         end
    end
    
    data = zeros(n_quad_nodes, n_hist_var+ hmsh.ndim);
    counter = 0;
    for ilev = 1:n_levels
         if (hmsh.nel_per_level(ilev) > 0)
            tot_nqn_lev =  hmsh.mesh_of_level(ilev).nqn * hmsh.nel_per_level(ilev);
            for idim = 1:hmsh.rdim                
                data(counter+1: counter+tot_nqn_lev, idim) = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), [tot_nqn_lev,1]);
            end
            data(counter+1: counter+tot_nqn_lev,idim+1:end) =  reshape (permute(eps_pl{ilev},[3,2,1]), [tot_nqn_lev,n_hist_var]);
            counter = counter+tot_nqn_lev;
         end         
    end

    % rhs = zeros (hspace.ndof, 1);
    % 
    % ndofs = 0;
    % for ilev = 1:hmsh.nlevels
    %   ndofs = ndofs + hspace.ndof_per_level(ilev);
    %   if (hmsh.nel_per_level(ilev) > 0)
    %     x = cell (hmsh.rdim, 1);
    %     for idim = 1:hmsh.rdim
    %       x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
    %     end
    %     sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
    %     b_lev = op_f_v (sp_lev, hmsh.msh_lev{ilev}, f(x{:}));
    % 
    %     dofs = 1:ndofs;
    %     rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
    %   end
    % end
    % 
    % n_quad_nodes = size(eps_pl,1) *  size(eps_pl,2);
    % data = zeros(n_quad_nodes, n_hist_var+ hmsh.ndim);
    % 
    % 
    % 
    % msh_val = msh_precompute(msh);
    % data(:,1:msh.ndim ) = reshape(msh_val.quad_nodes, msh.ndim, n_quad_nodes)';
    % data(:, msh.ndim+1:end) = reshape(permute(eps_pl, [3,2,1]),  n_hist_var , n_quad_nodes)';


    
    QI_coeff = get_QI_coeffs(hspace,hmsh,data);
    eps_pl_control_var(:,:) = QI_coeff;


end