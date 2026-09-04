function coeff = bezier_projection_levelwise (hspace, hmsh, eps_pl, measure)
% BEZIER_PROJECTION_LEVELWISE  Level-wise Bezier projection of quadrature-point history data onto a THB space
%
%   coeff = bezier_projection_levelwise (hspace, hmsh, eps_pl)
%   coeff = bezier_projection_levelwise (hspace, hmsh, eps_pl, measure)
%
% Each active function is fitted on its own level, where the element-local extraction is square,
% and the localized coefficients are averaged with the support weights of
% Thomas et al., CMAME 284 (2015), Eq. (61). Multi-level extraction: D'Angella et al., CMAME 328 (2018)

  % Tikhonov weight of the element-local fit. Keep at zero, any positive value destroys the convergence rate
  LAMBDA_LOC = 0;
  global BEZ_LAMBDA                                 % calibration override
  if (~isempty (BEZ_LAMBDA)); LAMBDA_LOC = BEZ_LAMBDA; end

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

  global BEZ_DIAG
  diag_on = ~isempty (BEZ_DIAG) && BEZ_DIAG;
  if (diag_on)
    d_cloc = 0; d_nel = 0; d_mincol = 0;
    d_ncontrib = zeros (hspace.ndof, 1);
    d_lev = zeros (hspace.ndof, 1);
  end

  % active functions are numbered level by level, coarsest first
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
    % w weights the element-local fit ('unit' gives the collocation of Hennig et al.), wphys the
    % support-weighted average, which must always use the physical measure (Thomas et al. Eq. (61))
    wphys = msh_lev.quad_weights .* msh_lev.jacdet;
    w = wphys;
    if (unit_measure); w = ones (msh_lev.nqn, msh_lev.nel); end
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

      % Weighted least-squares fit solved directly, never through the normal equations (they square cond(Ntp))
      sw = sqrt (wq);
      A  = sw .* Ntp;
      cn = vecnorm (A);
      if (max (cn) <= 0); continue; end             % zero-measure element
      good = cn > max (cn) * 1e-13;                 % functions that live here
      Ntp = Ntp(:, good); conn = conn(good); A = A(:, good);

      rhs = sw .* fq;
      if (LAMBDA_LOC > 0)
        nA = size (A, 2);
        scal = sqrt (LAMBDA_LOC) * norm (A, 'fro') / sqrt (nA);
        A = vertcat (A, scal * eye (nA));
        rhs = vertcat (rhs, zeros (nA, ncomp));
      end
      btp = A \ rhs;
      if (~all (isfinite (btp(:)))); continue; end

      % Average over the whole support where the extraction is invertible (Thomas et al. Lemmas 3.1
      % and 3.2), falling back to the own-level subset, which is square by construction
      CsubT_e  = CsubT(:, conn);
      rows_all = find (any (CsubT_e, 2));
      if (isempty (rows_all)); continue; end

      [rows, E] = extract_block (CsubT_e, rows_all);
      use_all = ~isempty (rows) && size (E, 2) <= size (E, 1) && cond (E) < 1e6;
      if (~use_all)
        [rows, E] = extract_block (CsubT_e, rows_all(rows_all >= lo & rows_all <= hi));
        if (isempty (rows)); continue; end
      end

      c_loc = E \ btp;

      if (diag_on)
        d_nel = d_nel + 1;
        d_mincol = max (d_mincol, max (abs (fq(:))));
        d_cloc = max (d_cloc, max (abs (btp(:))));
        d_ncontrib(rows) = d_ncontrib(rows) + 1;
        d_lev(rows) = ilev;
      end

      intN = (Ntp * E).' * wphys(:, iel);           % int_e T_A over this element
      num(rows, :) = num(rows, :) + intN .* c_loc;
      den(rows)    = den(rows)    + intN;
    end
  end

  coeff = zeros (hspace.ndof, ncomp);
  pos = den ~= 0;
  coeff(pos, :) = num(pos, :) ./ den(pos);

  if (diag_on)
    [cmax, iA] = max (max (abs (coeff), [], 2));
    fprintf (['[BEZ_DIAG] nlev=%d elems=%d  max|data|=%.3e  max|btp|=%.3e  ' ...
              'max|coeff|=%.3e @dof %d (lev %d, %d elems, den=%.3e, ' ...
              'den/max_den=%.2e)\n'], ...
             hmsh.nlevels, d_nel, d_mincol, d_cloc, cmax, iA, d_lev(iA), ...
             d_ncontrib(iA), den(iA), den(iA) / max (den));
  end
  if (any (~pos))
    error (['bezier_projection_levelwise: %d of %d functions received no data ' ...
            'from the active elements of their own level'], sum (~pos), hspace.ndof);
  end
end

function [rows, E] = extract_block (CsubT_e, rows_in)
% Extraction block of the given active functions on one element, truncated columns dropped
  if (isempty (rows_in)); rows = []; E = []; return; end
  E = full (CsubT_e(rows_in, :)).';
  keep = any (abs (E) > 0, 1);
  rows = rows_in(keep(:));
  E = E(:, keep);
end
