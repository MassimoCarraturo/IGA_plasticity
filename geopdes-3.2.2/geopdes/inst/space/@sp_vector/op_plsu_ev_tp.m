% OP_PLSU_EV_TP: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)), exploiting the tensor product structure.
%
%   mat = op_plsu_ev_tp (spu, spv, msh, lambda, mu);
%   [rows, cols, values] = op_plsu_ev_tp (spu, spv, msh, lambda, mu);
%
% INPUT:
%    
%   spu:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   lambda, mu: function handles to compute the Lame' coefficients
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2011, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [K, internal_energy, eps_pl, sigma] = op_plsu_ev_tp (space1, space2, msh, uhat, eps_pl, sigma, mu, kappa, sigma_y)

  for icomp = 1:space1.ncomp_param
    for idim = 1:msh.ndim
      size1 = size (space1.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      size2 = size (space2.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
        error ('One of the discrete spaces is not associated to the mesh')
      end
    end
  end

  K = spalloc (space2.ndof, space1.ndof, 5*space1.ndof); % tangent matrix
  internal_energy = zeros (space2.ndof, 1); % internal energy


  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    [A_temp, b_temp, eps_pl_temp, sigma_temp] = op_plsu_ev (sp1_col, sp2_col, msh_col, uhat, eps_pl(msh_col.elem_list,:,:), sigma(msh_col.elem_list,:,:), mu(x{:}), kappa(x{:}), sigma_y(x{:}));
    K = K + A_temp;
    internal_energy = internal_energy + b_temp;
    eps_pl(msh_col.elem_list,:,:) = eps_pl_temp;
    sigma(msh_col.elem_list,:,:) = sigma_temp;
  end
end