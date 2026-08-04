function coeff = bezier_projection_hier (hspace, hmsh, eps_pl, measure)
% BEZIER_PROJECTION_HIER  Bezier projection of quadrature-point history data
% onto a (truncated) hierarchical B-spline space.
%
%   coeff = bezier_projection_hier (hspace, hmsh, eps_pl)
%   coeff = bezier_projection_hier (hspace, hmsh, eps_pl, measure)
%
% INPUT
%   hspace : scalar hierarchical_space the data is projected onto
%   hmsh   : hierarchical_mesh carrying the quadrature data (must coincide
%            with the mesh underlying hspace; bisected projection meshes are
%            rejected by the caller)
%   eps_pl : cell array over levels, eps_pl{ilev} of size
%            [nel_per_level(ilev) x nqn x ncomp]
%   measure: (optional) 'physical' (default) or 'unit'.  'physical' weights
%            the element-local fit and the averaging by quad_weights.*jacdet,
%            i.e. an L2 projection on the physical element (Thomas et al.
%            Bezier projection).  'unit' drops the measure to identity, which
%            turns the element-local fit into the discrete (Euclidean) least
%            squares of Hennig et al. (CMAME 334, 2018, Eq. 41) with their
%            weighted local averaging (LLSQ_w, Eq. 33); this is the 'DLSQ_W'
%            projector.
%
% OUTPUT
%   coeff  : [hspace.ndof x ncomp] spline coefficients
%
% METHOD
%   Bezier projection in the sense of Thomas, Scott, Evans et al. (CMAME 284,
%   2015): an element-local L2 projection followed by a support-weighted
%   averaging of the localized coefficients,
%       c_A = sum_e w_A^e c_A^e,   w_A^e = int_e N_A / sum_e' int_e' N_A,
%   carried to the hierarchical setting through the multi-level (Bezier)
%   extraction operator of D'Angella, Kollmannsberger, Rank & Reali (CMAME
%   328, 2018): on an element of level l, every active THB function is a
%   linear combination of the level-l tensor-product B-splines (their
%   Secs. 3.1-3.6).  GeoPDEs stores exactly this operator level-wise in
%   hspace.Csub{l} (the global variant the paper discusses as ref. [23]);
%   restricting its rows to the tensor-product functions supported on the
%   element yields the local multi-level extraction operator.
%
%   The element-local projection is solved in the restricted THB basis
%   E_e * N_tp rather than in the Bernstein basis: the two span the same
%   local polynomial space, so the local L2 projection is identical (the
%   Bernstein detour of the references is a change of basis that cancels).
%   In the classical square case (exactly (p+1)^d functions on the element,
%   linearly independent) this reproduces Thomas et al. exactly; that case
%   is detected by a Cholesky factorization of the local Gramian and solved
%   directly.  On THB meshes an element can carry MORE functions than
%   dim Q_p -- the non-constant per-element DOF count that D'Angella et al.
%   note is inherent to the hierarchy (their Sec. 3.7) -- and the
%   restricted basis can then be linearly dependent.  In that case the
%   local system is solved with a truncated spectral pseudo-inverse: modes
%   of the Gramian with eigenvalue below RTOL*max(eig) are excluded from
%   the solve, and the SAME cutoff defines the local kernel, so no
%   near-singular mode is ever inverted.  The averaging is
%   dependence-aware: a function whose row norm over the kernel eigenspace
%   (a rotation-invariant measure) is significant is NOT identifiable from
%   that element -- its local coefficient is representation-dependent --
%   and receives zero averaging weight there.  Identifiable coefficients
%   are shared by every solution of the local system, so for a field that
%   lies in the spline space the averaging returns its exact coefficients.
%   A function identifiable on no element is an error, not a silent zero.
%
%   The local Gramian and load vector are assembled with the physical
%   measure (quad_weights .* jacdet), i.e. the projection is L2 on the
%   physical element, as in Thomas et al.

  if (nargin < 4 || isempty (measure))
    measure = 'physical';
  end
  if (~any (strcmpi (measure, {'physical', 'unit'})))
    error ('bezier_projection_hier: measure must be ''physical'' or ''unit''');
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
    error ('bezier_projection_hier: no active elements found');
  end

  % Spectral cutoff shared by the local solve and the identifiability test,
  % and rotation-invariant kernel-participation threshold.
  %
  % The cutoff must sit between two clusters: exact linear-dependence modes
  % of the restricted THB basis (s/smax at rounding level, <~ 1e-13) and
  % legitimate small-mass modes caused by near-degenerate geometry (the
  % eighth-sphere patch has a collapsed edge; functions squeezed by the map
  % produce s/smax down to ~1e-9 that carry real information).  1e-8 was
  % tried and misclassifies 54/462 functions on the sphere at max_level=3;
  % 1e-10 separates the clusters cleanly on every tested mesh (zero
  % non-identifiable functions across both benchmark studies).
  % The clusters above were measured at p = 2, where the local system carries
  % (p+1)^dim = 9 functions.  At p = 4 it carries 25 and the separation has to
  % be re-established, so the cutoff is overridable for calibration runs.
  tol_part = 1e-6;
  rtol_ker = 1e-10;
  ov = getenv ('BEZIER_RTOL_KER');
  if (~isempty (ov)); rtol_ker = str2double (ov); end

  num = zeros (hspace.ndof, ncomp);
  den = zeros (hspace.ndof, 1);

  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) == 0); continue; end

    msh_lev = hmsh.msh_lev{ilev};
    if (size (eps_pl{ilev}, 1) ~= msh_lev.nel || size (eps_pl{ilev}, 2) ~= msh_lev.nqn)
      error (['bezier_projection_hier: eps_pl{%d} is %dx%d but the mesh has ' ...
              '%d elements x %d quadrature points; data and scalar mesh do not coincide'], ...
             ilev, size (eps_pl{ilev}, 1), size (eps_pl{ilev}, 2), msh_lev.nel, msh_lev.nqn);
    end
    sp_lev  = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true);
    CsubT   = hspace.Csub{ilev}.';                  % last_dof(ilev) x ndof_TP(ilev): CSC column slicing
    if (unit_measure)
      w = ones (msh_lev.nqn, msh_lev.nel);          % Euclidean (discrete LSQ) measure
    else
      w = msh_lev.quad_weights .* msh_lev.jacdet;   % nqn x nel physical L2 measure
    end
    if (any (~isfinite (w(:))))
      error ('bezier_projection_hier: non-finite quadrature measure at level %d (degenerate geometry?)', ilev);
    end

    for iel = 1:msh_lev.nel
      conn = sp_lev.connectivity(:, iel);
      mask = conn > 0;
      conn = conn(mask);

      % local multi-level extraction: active THB functions on this element,
      % expressed on the level-l tensor-product functions supported on it
      CsubT_e = CsubT(:, conn);
      funcs   = find (any (CsubT_e, 2));            % global (cumulative) indices
      E       = full (CsubT_e(funcs, :)).';         % n_tp_loc x n_e

      Ntp  = sp_lev.shape_functions(:, mask, iel);  % nqn x n_tp_loc
      Nloc = Ntp * E;                               % nqn x n_e
      wq   = w(:, iel);

      G   = Nloc.' * (wq .* Nloc);                  % n_e x n_e local Gramian
      fq  = reshape (eps_pl{ilev}(iel, :, :), [], ncomp);
      rhs = Nloc.' * (wq .* fq);
      intN = Nloc.' * wq;                           % int_e N_A  (>= 0 for THB)

      % fast path: square extraction (the regular, non-interface case) gives a
      % genuine SPD B-spline Gramian; Cholesky both certifies and solves it
      n_e = numel (funcs);
      solved = false;
      if (n_e == nnz (mask))
        [R, flag] = chol ((G + G.') / 2);
        if (flag == 0)
          c_loc = R \ (R.' \ rhs);
          intN_id = intN;
          solved = true;
        end
      end

      if (~solved)
        % truncated spectral pseudo-inverse: one cutoff decides both which
        % modes are inverted and which functions are identifiable
        [V, S] = eig ((G + G.') / 2);
        s = diag (S);
        smax = max (s);
        if (~(smax > 0))                            % zero-measure element
          continue
        end
        keep = s > rtol_ker * smax;
        c_loc = V(:, keep) * ((V(:, keep).' * rhs) ./ s(keep));
        % rotation-invariant participation: row norm over the kernel eigenspace
        intN_id = intN;
        if (any (~keep))
          nonid = sqrt (sum (V(:, ~keep).^2, 2)) > tol_part;
          intN_id(nonid) = 0;
        end
      end

      num(funcs, :) = num(funcs, :) + intN_id .* c_loc;
      den(funcs)    = den(funcs)    + intN_id;
    end
  end

  coeff = zeros (hspace.ndof, ncomp);
  pos = den > 0;
  coeff(pos, :) = num(pos, :) ./ den(pos);

  if (any (~pos))
    error (['bezier_projection_hier: %d of %d functions identifiable on no ' ...
            'element of their support'], sum (~pos), hspace.ndof);
  end
end
