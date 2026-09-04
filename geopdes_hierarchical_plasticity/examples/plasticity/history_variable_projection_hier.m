function eps_pl_control_var = history_variable_projection_hier(hspace, hmsh, eps_pl, type_proj, hmsh_displ, lambda_override)
    % Project quadrature-point history variables onto the hierarchical scalar space
    % lambda_override forces the QI Tikhonov parameter (default 1e-9 on a bisected projection mesh, 0 otherwise)
    if nargin < 6
        lambda_override = [];
    end
    persistent localLS;
    if isempty(localLS)
        if exist('qi_local_ls_mex', 'file') == 3
            localLS = @getcoeff_localLS_Bspl_c;
        else
            localLS = @getcoeff_localLS_Bspl;
        end
    end

    if  strcmpi(type_proj, 'QI_C0')
        % the QI operator itself, the C^0 scalar space is selected upstream
        eps_pl_control_var = history_variable_projection_hier(hspace, hmsh, eps_pl, 'QI', hmsh_displ, lambda_override);

    else % QI or L2 projection
        % eps_pl lives on hmsh_displ, the projection itself uses the scalar space on hmsh
        hmsh_data = hmsh_displ;              % mesh that owns the data
        n_levels_data = hmsh_data.nlevels;
        n_hist_var = size(eps_pl{n_levels_data}, 3);
        eps_pl_control_var= zeros(hspace.ndof, n_hist_var);

        if  strcmpi(type_proj, 'L2')

            M = op_u_v_hier( hspace, hspace, hmsh);
            rhs = eval_rhs_l2_hier (hspace, hmsh, eps_pl, n_hist_var);
            % one factorisation for all components
            eps_pl_control_var = M \ rhs;

        elseif  strcmpi(type_proj, 'QI')

            n_quad_nodes = 0;
            for ilev = 1:n_levels_data
                 if (hmsh_data.nel_per_level(ilev) > 0)
                    n_quad_nodes = n_quad_nodes + hmsh_data.mesh_of_level(ilev).nqn * hmsh_data.nel_per_level(ilev);
                 end
            end

            % data(i,:) = point coordinates, f(i,:) = history components
            data = zeros(n_quad_nodes, hmsh_data.ndim);
            f = zeros(n_quad_nodes, n_hist_var);
            % Weight the local fits by the quadrature measure, otherwise the point density
            % imbalance across levels biases every function straddling two levels
            wquad = ones(n_quad_nodes, 1);
            counter = 0;
            for ilev = 1:n_levels_data
                 if (hmsh_data.nel_per_level(ilev) > 0)
                    tot_nqn_lev =  hmsh_data.mesh_of_level(ilev).nqn * hmsh_data.nel_per_level(ilev);

                    quad_nodes = my_msh_evaluate_qn (hmsh_data.mesh_of_level(ilev), hmsh_data.active{ilev});
                    for idim = 1:hmsh_data.rdim
                        data(counter+1: counter+tot_nqn_lev, idim) =  reshape (quad_nodes(idim,:,:), [tot_nqn_lev,1]);
                    end
                    f(counter+1: counter+tot_nqn_lev,:) =  reshape (permute(eps_pl{ilev},[3,2,1]), [n_hist_var, tot_nqn_lev])';
                    mlev = msh_evaluate_element_list (hmsh_data.mesh_of_level(ilev), hmsh_data.active{ilev});
                    wlev = mlev.quad_weights .* mlev.jacdet;       % nqn x nel
                    wquad(counter+1: counter+tot_nqn_lev) = reshape (wlev, [tot_nqn_lev, 1]);
                            counter = counter+tot_nqn_lev;
                 end
            end

            % Tikhonov term only when the projection mesh is finer than the data mesh.
            % The penalty scales like h^-4, so a fixed lambda on coincident meshes destroys convergence
            if isempty(lambda_override)
                lambda_qi = (hmsh.nlevels > hmsh_displ.nlevels) * 1e-9;
            else
                lambda_qi = lambda_override;
            end
            % Speleers and Manni, Numer. Math. 132 (2016), Theorem 4: local LS is not a
            % projector on a hierarchy, so the Speleers QI is used when nlevels > 1
            if (hmsh.nlevels > 1)
                QI_coeff = getcoeff_qi_speleers(hspace, hmsh, eps_pl);
            else
                QI_coeff = localLS(hspace,hmsh,data,f,lambda_qi,wquad.');
            end

            eps_pl_control_var(:,:) = QI_coeff;

        elseif  any (strcmpi(type_proj, {'BEZIER', 'BEZIER_C0', 'DLSQ', 'DLSQ_W'}))

            % eps_pl rows must follow the scalar mesh element list
            if (hmsh.nlevels ~= hmsh_displ.nlevels || ...
                ~isequal (hmsh.nel_per_level, hmsh_displ.nel_per_level))
                error (['history_variable_projection_hier: %s requires the scalar mesh ' ...
                        'to coincide with the displacement mesh (no bisected projection mesh)'], upper(type_proj));
            end
            switch upper(type_proj)
                case 'DLSQ'     % Hennig et al. 2018, Eqs. (41)-(42)
                    eps_pl_control_var = hennig_dlsq_projection_hier (hspace, hmsh, eps_pl);
                case 'DLSQ_W'   % Hennig et al. 2018, Eq. (33): Bezier projector with the unit measure
                    eps_pl_control_var = bezier_projection_levelwise (hspace, hmsh, eps_pl, 'unit');
                otherwise       % Bezier projection via hspace.Csub
                    eps_pl_control_var = bezier_projection_levelwise (hspace, hmsh, eps_pl);
            end

        else
            error ('history_variable_projection_hier: unknown type_proj ''%s''', type_proj);
        end
    end
end


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


