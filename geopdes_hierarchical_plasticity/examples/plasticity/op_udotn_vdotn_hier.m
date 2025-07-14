function varargout = op_udotn_vdotn_hier (hspace, hmsh, coeff, iside)



  K = spalloc (hspace.ndof, hspace.ndof, 3*hspace.ndof);



  ndofs = 0;
  for ilev = 1:hmsh.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (ilev <= numel(hmsh.boundary(iside).nel_per_level) && hmsh.boundary(iside).nel_per_level(ilev) > 0)
      x = cell (hmsh.rdim, 1);

      hmsh.mesh_of_level(ilev).boundary = hmsh.boundary(1).mesh_of_level(ilev);
      for is = 2:numel(hmsh.boundary)           
          hmsh.mesh_of_level(ilev).boundary(is) = hmsh.boundary(iside).mesh_of_level(ilev)  ;
      end

      msh_side = msh_eval_boundary_side (hmsh.mesh_of_level(ilev), iside, hmsh.boundary(iside).active{ilev});
      msh_side_from_interior = msh_boundary_side_from_interior (hmsh.mesh_of_level(ilev), iside);

      sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
      msh_side_from_interior_struct = msh_evaluate_element_list(msh_side_from_interior, hmsh.boundary(iside).active{ilev});
      sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);

      for idim = 1:hmsh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end


      K_lev = op_udotn_vdotn (sp_bnd, sp_bnd,msh_side, coeff(x{:}));
      dofs = 1:ndofs;
      K(dofs,dofs) = K(dofs,dofs) + hspace.Csub{ilev}.' * K_lev * hspace.Csub{ilev};
    end
  end

  if (nargout == 1)
    varargout{1} = K;
  elseif (nargout == 3)
    [rows, cols, vals] = find (K);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end


% rhs = zeros (hspace.ndof, 1);
% 
%   ndofs = 0;
%   for ilev = 1:hmsh.nlevels
%     ndofs = ndofs + hspace.ndof_per_level(ilev);
%     if (ilev <= numel(hmsh.boundary(iside).nel_per_level) && hmsh.boundary(iside).nel_per_level(ilev) > 0)
%       x = cell (hmsh.rdim, 1);
%       %%% use this trick to store the boundary information on the fly
%       %%% inside mesh_of_level
%       hmsh.mesh_of_level(ilev).boundary = hmsh.boundary(iside).mesh_of_level(ilev);
%       msh_side = msh_eval_boundary_side (hmsh.mesh_of_level(ilev), iside, hmsh.boundary(iside).active{ilev});
%       msh_side_from_interior = msh_boundary_side_from_interior (hmsh.mesh_of_level(ilev), iside);
% 
%       sp_bnd = hspace.space_of_level(ilev).constructor (msh_side_from_interior);
%       msh_side_from_interior_struct = msh_evaluate_element_list(msh_side_from_interior, hmsh.boundary(iside).active{ilev});
%       sp_bnd = sp_evaluate_element_list (sp_bnd, msh_side_from_interior_struct, 'value', true, 'gradient', true);
% 
%       for idim = 1:hmsh.rdim
%         x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
%       end
% 
%       % b_lev = op_gradv_n_f (sp_bnd, msh_side, f(x{:}, iside));
%       b_lev = op_pn_v (sp_bnd, msh_side, f(x{:}) );
% 
%       dofs = 1:ndofs;
%       rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
%     end
%   end
% 
% end