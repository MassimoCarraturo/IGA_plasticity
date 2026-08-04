function eps_pl_control_var = history_variable_projection_hier(hspace, hmsh, eps_pl, type_proj, hmsh_displ, lambda_override)
    % libqi backend (C/OpenMP) for the local-LS B-spline projection.
    % Falls back to the MATLAB reference if the MEX is unavailable.
    %
    % lambda_override (optional): forces the Tikhonov parameter of the QI
    % local least-squares systems.  When omitted, lambda is chosen
    % automatically (1e-9 on a finer projection mesh, 0 otherwise).  Used to
    % reproduce the unregularised QI-h runs of the article.
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
        n_hist_var = size(eps_pl{end}, 3);
        eps_pl_control_var= zeros(hspace.ndof, n_hist_var);

        n_quad_nodes = 0;
        for ilev = 1:hmsh_displ.nlevels
             if (hmsh_displ.nel_per_level(ilev) > 0)
                n_quad_nodes = n_quad_nodes + hmsh_displ.mesh_of_level(ilev).nqn * hmsh_displ.nel_per_level(ilev);
             end
        end

        data = zeros(n_quad_nodes, hmsh_displ.ndim);
        f = zeros(n_quad_nodes, n_hist_var);
        counter = 0;
        for ilev = 1:hmsh_displ.nlevels
             if (hmsh_displ.nel_per_level(ilev) > 0)
                tot_nqn_lev =  hmsh_displ.mesh_of_level(ilev).nqn * hmsh_displ.nel_per_level(ilev);

                % msh_lev = msh_evaluate_element_list (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
                quad_nodes = my_msh_evaluate_qn (hmsh_displ.mesh_of_level(ilev), hmsh_displ.active{ilev});
                for idim = 1:hmsh_displ.rdim                
                    % data(counter+1: counter+tot_nqn_lev, idim) = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), [tot_nqn_lev,1]);                  
                    data(counter+1: counter+tot_nqn_lev, idim) =  reshape (quad_nodes(idim,:,:), [tot_nqn_lev,1]);
                end
                f(counter+1: counter+tot_nqn_lev,:) =  reshape (permute(eps_pl{ilev},[3,2,1]), [n_hist_var, tot_nqn_lev])';
                counter = counter+tot_nqn_lev;
             end         
        end

        % compute spline coefficients
        QI_coeff = localLS(hspace,hmsh,data,f,1e-9);

        % store in output
        eps_pl_control_var(:,:) = QI_coeff;

    
    else % QI or L2 projection
        % When num_bisections > 0, hmsh (scalar mesh) may have more levels
        % than hmsh_displ (primal mesh).  The quadrature data (eps_pl) lives
        % on hmsh_displ, so we use hmsh_displ for the data assembly loop.
        % The scalar space/mesh is used for the projection itself.
        hmsh_data = hmsh_displ;              % mesh that owns the data
        n_levels_data = hmsh_data.nlevels;
        n_hist_var = size(eps_pl{n_levels_data}, 3);
        eps_pl_control_var= zeros(hspace.ndof, n_hist_var);

        if  strcmpi(type_proj, 'L2')

            M = op_u_v_hier( hspace, hspace, hmsh);
            rhs = eval_rhs_l2_hier (hspace, hmsh, eps_pl, n_hist_var);
            % one factorisation for all components; column-by-column M\rhs(:,i)
            % refactorises M every time
            eps_pl_control_var = M \ rhs;

        elseif  strcmpi(type_proj, 'QI')

            % compute number of quadrature points (on the data mesh)
            n_quad_nodes = 0;
            for ilev = 1:n_levels_data
                 if (hmsh_data.nel_per_level(ilev) > 0)
                    n_quad_nodes = n_quad_nodes + hmsh_data.mesh_of_level(ilev).nqn * hmsh_data.nel_per_level(ilev);
                 end
            end

            % organize point coordinates in data(i,:) = (coord_x, coord_y, ...)
            % organize plastic variables in f(i,:) = (plastic_var_1, plastic_var2, ...)
            data = zeros(n_quad_nodes, hmsh_data.ndim);
            f = zeros(n_quad_nodes, n_hist_var);
            counter = 0;
            for ilev = 1:n_levels_data
                 if (hmsh_data.nel_per_level(ilev) > 0)
                    tot_nqn_lev =  hmsh_data.mesh_of_level(ilev).nqn * hmsh_data.nel_per_level(ilev);

                    quad_nodes = my_msh_evaluate_qn (hmsh_data.mesh_of_level(ilev), hmsh_data.active{ilev});
                    for idim = 1:hmsh_data.rdim
                        data(counter+1: counter+tot_nqn_lev, idim) =  reshape (quad_nodes(idim,:,:), [tot_nqn_lev,1]);
                    end
                    f(counter+1: counter+tot_nqn_lev,:) =  reshape (permute(eps_pl{ilev},[3,2,1]), [n_hist_var, tot_nqn_lev])';
                    counter = counter+tot_nqn_lev;
                 end
            end

            % compute spline coefficients (project onto scalar space).
            % Regularize the per-DOF local LS only when the projection
            % (scalar) mesh is finer than the data mesh — i.e. num_bisections>0,
            % detected via extra hierarchical levels.  There the local systems
            % are under-determined (sparse/one-sided support coverage), so a
            % small Tikhonov term (matching the QI_C0 path) stabilizes them.
            % For num_bisections==0 the meshes coincide and lambda stays 0,
            % leaving plain QI / QI_p3 / QI_graded results unchanged.
            if isempty(lambda_override)
                lambda_qi = (hmsh.nlevels > hmsh_displ.nlevels) * 1e-9;
            else
                lambda_qi = lambda_override;
            end
            QI_coeff = localLS(hspace,hmsh,data,f,lambda_qi);

            % store in output
            eps_pl_control_var(:,:) = QI_coeff;

        elseif  any (strcmpi(type_proj, {'BEZIER', 'BEZIER_C0'}))

            % Bezier projection (element-local L2 + support-weighted averaging,
            % via the multi-level extraction operator stored in hspace.Csub).
            % BEZIER_C0 is the identical operator applied to a C^0 scalar space
            % (interior knots of multiplicity p); the space is selected upstream
            % in adaptivity_J2_plasticity, so nothing changes here.
            % The data must live on the same mesh the scalar space is built on:
            % a bisected projection mesh (num_bisections > 0) is not supported,
            % and the element ordering of the two meshes must coincide because
            % eps_pl rows are indexed by the scalar mesh's element list.
            if (hmsh.nlevels ~= hmsh_displ.nlevels || ...
                ~isequal (hmsh.nel_per_level, hmsh_displ.nel_per_level))
                error (['history_variable_projection_hier: BEZIER requires the scalar mesh ' ...
                        'to coincide with the displacement mesh (no bisected projection mesh)']);
            end
            eps_pl_control_var = bezier_projection_levelwise (hspace, hmsh, eps_pl);

        elseif  strcmpi(type_proj, 'DLSQ')

            % Global discrete least-squares projection (Hennig et al. 2018,
            % Eqs. 41-42): quadrature-weight-free Euclidean fit at the QPs.
            % Same mesh-coincidence requirement as BEZIER.
            if (hmsh.nlevels ~= hmsh_displ.nlevels || ...
                ~isequal (hmsh.nel_per_level, hmsh_displ.nel_per_level))
                error (['history_variable_projection_hier: DLSQ requires the scalar mesh ' ...
                        'to coincide with the displacement mesh (no bisected projection mesh)']);
            end
            eps_pl_control_var = hennig_dlsq_projection_hier (hspace, hmsh, eps_pl);

        elseif  strcmpi(type_proj, 'DLSQ_W')

            % Element-local discrete least-squares with weighted averaging
            % (Hennig et al. 2018, LLSQ_w, Eq. 33): the Bezier projector with
            % the identity (unit) measure instead of the physical L2 measure.
            if (hmsh.nlevels ~= hmsh_displ.nlevels || ...
                ~isequal (hmsh.nel_per_level, hmsh_displ.nel_per_level))
                error (['history_variable_projection_hier: DLSQ_W requires the scalar mesh ' ...
                        'to coincide with the displacement mesh (no bisected projection mesh)']);
            end
            eps_pl_control_var = bezier_projection_levelwise (hspace, hmsh, eps_pl, 'unit');

        else
            error ('history_variable_projection_hier: unknown type_proj ''%s''', type_proj);
        end
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


