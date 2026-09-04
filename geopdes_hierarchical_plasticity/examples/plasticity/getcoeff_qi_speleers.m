% GETCOEFF_QI_SPELEERS  Hierarchical quasi-interpolant projector onto the THB-spline space
%
%   QI_coeff = getcoeff_qi_speleers (hspace, hmsh, eps_pl)
%
% eps_pl: cell per level, (nel x nqn x ncomp) over the active elements of that level
% Speleers and Manni, Numer. Math. 132 (2016), Eq. (22) and Theorem 4
% Per-cell candidates are averaged with the weights of Thomas et al., CMAME 284 (2015), Eq. (61)
% A single cell per coefficient still reproduces the space but blows up at the plastic front

function QI_coeff = getcoeff_qi_speleers (hspace, hmsh, eps_pl)

  ncomp = 0;
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) > 0)
      ncomp = size (eps_pl{ilev}, 3);
      break
    end
  end
  if (ncomp == 0)
    error ('getcoeff_qi_speleers: no active elements found');
  end

  QI_coeff = zeros (hspace.ndof, ncomp);
  done     = false (hspace.ndof, 1);

  last_of_lev  = cumsum (hspace.ndof_per_level(:).');
  first_of_lev = [1, last_of_lev(1:end-1) + 1];

  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) == 0);        continue; end
    if (hspace.ndof_per_level(ilev) == 0);     continue; end

    msh_lev = hmsh.msh_lev{ilev};
    sp_lev  = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true);

    % gdof is zero where the level-l function is not active
    gdof = zeros (hspace.space_of_level(ilev).ndof, 1);
    act  = hspace.active{ilev}(:);
    gdof(act) = first_of_lev(ilev) : last_of_lev(ilev);

    acc  = zeros (hspace.ndof, ncomp);
    wsum = zeros (hspace.ndof, 1);
    wphys = msh_lev.quad_weights .* msh_lev.jacdet;

    for iel = 1:msh_lev.nel
      conn = sp_lev.connectivity(:, iel);
      mask = conn > 0;
      if (~any (mask)); continue; end
      cn = conn(mask);

      g  = gdof(cn);
      ok = g > 0;
      if (~any (ok)); continue; end

      % Speleers and Manni Eqs. (38) and (40), least squares when nquad > p+1
      N  = sp_lev.shape_functions(:, mask, iel);
      fq = reshape (eps_pl{ilev}(iel, :, :), [], ncomp);
      if (size (N, 1) == size (N, 2) && rcond (N) < eps); continue; end
      c = N \ fq;
      if (~all (isfinite (c(:)))); continue; end

      wA = N.' * wphys(:, iel);

      gi = g(ok);  wi = max (wA(ok), 0);
      acc(gi, :) = acc(gi, :) + wi .* c(ok, :);
      wsum(gi)   = wsum(gi) + wi;
    end

    have = wsum > 0;
    QI_coeff(have, :) = acc(have, :) ./ wsum(have);
    done(have) = true;
  end

  % Functions missed above (Eq. (9) can fail after admissible coarsening) are
  % served through Csub from any active element that sees them, outside Theorem 4
  if (~any (done)); return; end
  for ilev = 1:hmsh.nlevels
    if (all (done)); break; end
    if (hmsh.nel_per_level(ilev) == 0); continue; end
    msh_lev = hmsh.msh_lev{ilev};
    sp_lev  = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true);
    CsubT   = hspace.Csub{ilev}.';
    for iel = 1:msh_lev.nel
      if (all (done)); break; end
      conn = sp_lev.connectivity(:, iel);
      mask = conn > 0;
      if (~any (mask)); continue; end
      cn = conn(mask);
      CsubT_e = CsubT(:, cn);
      rows = find (any (CsubT_e, 2));
      rows = rows(~done(rows));
      if (isempty (rows)); continue; end

      N  = sp_lev.shape_functions(:, mask, iel);
      fq = reshape (eps_pl{ilev}(iel, :, :), [], ncomp);
      E  = full (CsubT_e(rows, :)).';
      keep = any (abs (E) > 0, 1);
      rows = rows(keep(:));  E = E(:, keep);
      if (isempty (rows) || size (E, 2) > size (E, 1)); continue; end
      btp = N \ fq;
      c   = E \ btp;
      if (~all (isfinite (c(:)))); continue; end
      QI_coeff(rows, :) = c;
      done(rows) = true;
    end
  end

  if (~all (done))
    warning ('getcoeff_qi_speleers:unassigned', ...
             '%d of %d active functions received no value and stay at zero', ...
             sum (~done), hspace.ndof);
  end
end
