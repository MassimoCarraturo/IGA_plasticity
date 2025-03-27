function eps_pl_control_var = history_variable_projection_hier(hspace, hmsh, eps_pl, type_proj)

    n_levels = hmsh.nlevels;
    n_hist_var = size(eps_pl{n_levels}, 3);
    eps_pl_control_var= zeros(hspace.ndof, n_hist_var);
    
    if  strcmpi(type_proj, 'L2')
        
        M = op_u_v_hier( hspace, hspace, hmsh);  
        rhs = eval_rhs_l2_hier (hspace, hmsh, eps_pl, n_hist_var);
        for i=1:n_hist_var           
            eps_pl_control_var(:,i) = M\rhs(:,i);        
        end
    
    elseif  strcmpi(type_proj, 'QI')
        if hspace.ncomp == 3
            print('QI for 2D only!!!!!!!')
            return
        end
        
        % compute number of quadrature points
        n_quad_nodes = 0;
        for ilev = 1:n_levels
             if (hmsh.nel_per_level(ilev) > 0)
                n_quad_nodes = n_quad_nodes + hmsh.mesh_of_level(ilev).nqn * hmsh.nel_per_level(ilev);
             end
        end
        
        % organize plastic variables in data(i,:) = (coord_x, coord_y, platic_var_1, platic_var2, ...)
        data = zeros(n_quad_nodes, n_hist_var+ hmsh.ndim);
        counter = 0;
        for ilev = 1:n_levels
             if (hmsh.nel_per_level(ilev) > 0)
                tot_nqn_lev =  hmsh.mesh_of_level(ilev).nqn * hmsh.nel_per_level(ilev);

                % msh_lev = msh_evaluate_element_list (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
                quad_nodes = my_msh_evaluate_qn (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
                for idim = 1:hmsh.rdim                
                    % data(counter+1: counter+tot_nqn_lev, idim) = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), [tot_nqn_lev,1]);                  
                    data(counter+1: counter+tot_nqn_lev, idim) =  reshape (quad_nodes(idim,:,:), [tot_nqn_lev,1]);
                end
                data(counter+1: counter+tot_nqn_lev,idim+1:end) =  reshape (permute(eps_pl{ilev},[3,2,1]), [n_hist_var, tot_nqn_lev])';
                counter = counter+tot_nqn_lev;
             end         
        end
    
        % compute spline coefficients
        QI_coeff = get_QI_coeffs(hspace,hmsh,data);
        % store in output
        eps_pl_control_var(:,:) = QI_coeff;
    
    
    end
end

%---------------------------------------------------------------------------------------------------
% subroutines for L2 projection

function rhs = eval_rhs_l2_hier (hspv, hmsh, eps_pl,n_hist_var)

    rhs = zeros (hspv.ndof, n_hist_var);
    ndofs_v = 0;
    for ilev = 1:hmsh.nlevels
        ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
        if (hmsh.nel_per_level(ilev) > 0)
    
            spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), hmsh.msh_lev{ilev});
            dofs_v = 1:ndofs_v;     
            
            for ivar = 1:n_hist_var
                rhs_lev = op_f_v (spv_lev, hmsh.msh_lev{ilev}, eps_pl{ilev}(:,:,ivar)');         
                rhs(dofs_v,ivar) = rhs(dofs_v,ivar) + hspv.Csub{ilev}.'*rhs_lev;
            end            
        end
    end
end


