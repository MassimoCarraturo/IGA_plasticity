function coeff = hennig_dlsq_projection_hier (hspace, hmsh, eps_pl)
% HENNIG_DLSQ_PROJECTION_HIER  Global discrete least-squares projection of
% quadrature-point history data onto a (truncated) hierarchical B-spline space.
%
%   coeff = hennig_dlsq_projection_hier (hspace, hmsh, eps_pl)
%
% INPUT
%   hspace : scalar hierarchical_space the data is projected onto
%   hmsh   : hierarchical_mesh carrying the quadrature data (must coincide
%            with the mesh underlying hspace; bisected projection meshes are
%            rejected by the caller)
%   eps_pl : cell array over levels, eps_pl{ilev} of size
%            [nel_per_level(ilev) x nqn x ncomp]
%
% OUTPUT
%   coeff  : [hspace.ndof x ncomp] spline coefficients
%
% METHOD
%   The global discrete least-squares (LSQ) projector of Hennig, Ambati,
%   De Lorenzis & Kaestner (CMAME 334, 2018, Sec. 3.3.1(ii), Eqs. 41-42):
%   the control values H[A] minimise the EUCLIDEAN residual at the
%   quadrature points,
%       min_H[A]  || H[Q] - T H[A] ||_2^2 ,   T_IJ = N_J(x_I),
%   where N_J are the hierarchical basis functions and x_I the physical
%   quadrature-point coordinates.  Unlike the global L2 projection this is
%   quadrature-weight-free (identity measure), and unlike the local Bezier
%   projector it solves one global system.  The collocation matrix T is
%   assembled level-wise through the multi-level extraction operator stored
%   in hspace.Csub, then the overdetermined system is solved by sparse QR
%   least squares (coeff = T \ H), which is the exact minimiser of Eq. 41.

  % ---- detect number of history components ----
  ncomp = 0;
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      ncomp = size (eps_pl{ilev}, 3);
      break
    end
  end
  if (ncomp == 0)
    error ('hennig_dlsq_projection_hier: no active elements found');
  end

  % ---- total number of quadrature points across active elements ----
  ntot = 0;
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      ntot = ntot + hmsh.msh_lev{ilev}.nel * hmsh.msh_lev{ilev}.nqn;
    end
  end

  % ---- assemble collocation matrix T (ntot x ndof) and data H (ntot x ncomp) ----
  % Triplet buffers for T, filled per element and consistent with the row
  % ordering used for H.
  Ii = cell (1, hmsh.nlevels);
  Jj = cell (1, hmsh.nlevels);
  Vv = cell (1, hmsh.nlevels);
  H  = zeros (ntot, ncomp);
  last_dof = cumsum (hspace.ndof_per_level);
  rowoff = 0;

  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) == 0); continue; end

    msh_lev = hmsh.msh_lev{ilev};
    if (size (eps_pl{ilev}, 1) ~= msh_lev.nel || size (eps_pl{ilev}, 2) ~= msh_lev.nqn)
      error (['hennig_dlsq_projection_hier: eps_pl{%d} is %dx%d but the mesh has ' ...
              '%d elements x %d quadrature points; data and scalar mesh do not coincide'], ...
             ilev, size (eps_pl{ilev}, 1), size (eps_pl{ilev}, 2), msh_lev.nel, msh_lev.nqn);
    end
    sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true);
    CsubT  = hspace.Csub{ilev}.';                   % last_dof(ilev) x ndof_TP(ilev)
    nqn    = msh_lev.nqn;

    il = cell (1, msh_lev.nel);
    jl = cell (1, msh_lev.nel);
    vl = cell (1, msh_lev.nel);

    for iel = 1:msh_lev.nel
      conn = sp_lev.connectivity(:, iel);
      mask = conn > 0;
      conn = conn(mask);

      % active THB functions on this element and their local extraction
      CsubT_e = CsubT(:, conn);
      funcs   = find (any (CsubT_e, 2));            % global (cumulative) dof indices
      E       = full (CsubT_e(funcs, :)).';         % n_tp_loc x n_e
      Ntp     = sp_lev.shape_functions(:, mask, iel);   % nqn x n_tp_loc
      Nloc    = Ntp * E;                            % nqn x n_e : THB values at QPs

      rows = rowoff + (1:nqn).';
      ne   = numel (funcs);
      il{iel} = repmat (rows, ne, 1);
      jl{iel} = reshape (repmat (funcs(:).', nqn, 1), [], 1);
      vl{iel} = Nloc(:);

      H(rows, :) = reshape (eps_pl{ilev}(iel, :, :), nqn, ncomp);
      rowoff = rowoff + nqn;
    end

    Ii{ilev} = vertcat (il{:});
    Jj{ilev} = vertcat (jl{:});
    Vv{ilev} = vertcat (vl{:});
  end

  T = sparse (vertcat (Ii{:}), vertcat (Jj{:}), vertcat (Vv{:}), ntot, hspace.ndof);

  % ---- solve the discrete least-squares system (Eq. 41) ----
  % Sparse mldivide on an overdetermined system performs a QR least-squares
  % solve, i.e. exactly argmin ||T coeff - H||_2.
  coeff = T \ H;
end
