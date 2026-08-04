function coeff = bezier_projection_levelwise (hspace, hmsh, eps_pl, measure)
% BEZIER_PROJECTION_LEVELWISE  Level-wise Bezier projection of quadrature-point
% history data onto a (truncated) hierarchical B-spline space.
%
%   coeff = bezier_projection_levelwise (hspace, hmsh, eps_pl)
%   coeff = bezier_projection_levelwise (hspace, hmsh, eps_pl, measure)
%
% Arguments are those of BEZIER_PROJECTION_HIER.
%
% WHY THIS EXISTS
%   BEZIER_PROJECTION_HIER fits, on each element, all active THB functions
%   overlapping it.  Restricted to one element every such function is a
%   polynomial of degree <= p per direction, so at most (p+1)^dim of them are
%   independent; the hierarchy routinely puts more than that on an element, and
%   the count grows with the degree and with the number of levels meeting there.
%   The element-local system is then singular by construction.  At p = 2 enough
%   functions still survive the identifiability test for the support-weighted
%   average to reconstruct them all, but from p = 3 upward some functions are
%   identifiable on no element of their support and the operator has no value to
%   assign them.
%
% THE FIX
%   Fit each function on its OWN level.  On an active element of level l the
%   level-l tensor-product functions supported on it are exactly (p+1)^dim and
%   linearly independent, so the element-local L2 projection is a square SPD
%   solve at any degree.  This is precisely the classical Bezier projection of
%   Thomas et al. (CMAME 284, 2015), applied on the level where extraction is
%   square, with the multi-level extraction of D'Angella et al. (CMAME 328,
%   2018) used only to express the active THB functions of that level in the
%   level-l tensor-product basis.
%
%   The scheme is well posed because of the THB activation rule itself: a
%   function of level l is active only if its support lies in Omega^l but NOT in
%   Omega^{l+1}, so part of its support is always covered by active level-l
%   elements.  Every active function therefore receives data from its own level,
%   and no function can be left unassigned.
%
%   Localized coefficients are averaged with the same support weights as the
%   element-local operator,
%       c_A = sum_e w_A^e c_A^e,   w_A^e = int_e T_A / sum_e' int_e' T_A,
%   the sums running over the active elements of level(A) only.

  if (nargin < 4 || isempty (measure))
    measure = 'physical';
  end
  if (~any (strcmpi (measure, {'physical', 'unit'})))
    error ('bezier_projection_levelwise: measure must be ''physical'' or ''unit''');
  end
  unit_measure = strcmpi (measure, 'unit');

  ncomp = 0;
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      ncomp = size (eps_pl{ilev}, 3);
      break
    end
  end
  if (ncomp == 0)
    error ('bezier_projection_levelwise: no active elements found');
  end

  num = zeros (hspace.ndof, ncomp);
  den = zeros (hspace.ndof, 1);

  % global index range of the active functions of each level (they are numbered
  % level by level, coarsest first)
  last_of_lev  = cumsum (hspace.ndof_per_level(:).');
  first_of_lev = [1, last_of_lev(1:end-1) + 1];

  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) == 0); continue; end
    if (hspace.ndof_per_level(ilev) == 0); continue; end

    msh_lev = hmsh.msh_lev{ilev};
    if (size (eps_pl{ilev}, 1) ~= msh_lev.nel || size (eps_pl{ilev}, 2) ~= msh_lev.nqn)
      error (['bezier_projection_levelwise: eps_pl{%d} is %dx%d but the mesh has ' ...
              '%d elements x %d quadrature points'], ilev, ...
             size (eps_pl{ilev}, 1), size (eps_pl{ilev}, 2), msh_lev.nel, msh_lev.nqn);
    end

    sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true);
    CsubT  = hspace.Csub{ilev}.';                   % active THB x level-l TP
    if (unit_measure)
      w = ones (msh_lev.nqn, msh_lev.nel);
    else
      w = msh_lev.quad_weights .* msh_lev.jacdet;
    end
    if (any (~isfinite (w(:))))
      error ('bezier_projection_levelwise: non-finite quadrature measure at level %d', ilev);
    end

    lo = first_of_lev(ilev); hi = last_of_lev(ilev);

    for iel = 1:msh_lev.nel
      conn = sp_lev.connectivity(:, iel);
      mask = conn > 0;
      conn = conn(mask);

      Ntp = sp_lev.shape_functions(:, mask, iel);   % nqn x (p+1)^dim, independent
      wq  = w(:, iel);
      fq  = reshape (eps_pl{ilev}(iel, :, :), [], ncomp);

      % square, SPD element-local projection on the level-l tensor-product basis
      Gtp = Ntp.' * (wq .* Ntp);
      [Rc, flag] = chol ((Gtp + Gtp.') / 2);
      if (flag ~= 0); continue; end                 % zero-measure element
      btp = Rc \ (Rc.' \ (Ntp.' * (wq .* fq)));     % localized TP coefficients

      % active THB functions OF THIS LEVEL supported on this element, expressed
      % in the level-l tensor-product basis
      CsubT_e = CsubT(:, conn);
      rows    = find (any (CsubT_e, 2));
      rows    = rows(rows >= lo & rows <= hi);
      if (isempty (rows)); continue; end
      E = full (CsubT_e(rows, :)).';                % (p+1)^dim x n_lev

      keep = any (abs (E) > 0, 1);                  % truncated away on this element
      rows = rows(keep(:));
      E    = E(:, keep);
      if (isempty (rows)); continue; end

      % E has full column rank (distinct truncated B-splines of one level are
      % independent on an element), so this is a well-posed least-squares solve
      c_loc = E \ btp;

      intN = (Ntp * E).' * wq;                      % int_e T_A over this element
      num(rows, :) = num(rows, :) + intN .* c_loc;
      den(rows)    = den(rows)    + intN;
    end
  end

  coeff = zeros (hspace.ndof, ncomp);
  pos = den ~= 0;
  coeff(pos, :) = num(pos, :) ./ den(pos);
  if (any (~pos))
    error (['bezier_projection_levelwise: %d of %d functions received no data ' ...
            'from the active elements of their own level'], sum (~pos), hspace.ndof);
  end
end
