function eps_new = cpt_transfer_hier (hmsh_old, eps_old, hmsh_new)
% CPT_TRANSFER_HIER  Closest-point transfer of quadrature-point history data
% between two hierarchical meshes of the same geometry.
%
%   eps_new = cpt_transfer_hier (hmsh_old, eps_old, hmsh_new)
%
% INPUT
%   hmsh_old : hierarchical_mesh the data currently lives on
%   eps_old  : cell over old levels, eps_old{ilev} = [nel x nqn x ncomp]
%   hmsh_new : hierarchical_mesh to transfer the data onto
%
% OUTPUT
%   eps_new  : cell over new levels, eps_new{ilev} = [nel x nqn x ncomp]
%
% METHOD
%   Closest-point transfer (CPT) of Hennig, Ambati, De Lorenzis & Kaestner
%   (CMAME 334, 2018, Sec. 3.3.1(i)): each new quadrature point takes the
%   value of the geometrically nearest old quadrature point.  Hennig restrict
%   the search to the parent element for efficiency; for nested refine/coarsen
%   the global nearest point lies in that parent anyway, so we search over all
%   old quadrature points by physical coordinate.  This is the simplest and
%   least accurate transfer: exact only for element-wise constant data, and at
%   most first order otherwise.

  % ---- gather all old quadrature points: physical coords + values ----
  ncomp = 0;
  ntot_old = 0;
  for ilev = 1:hmsh_old.nlevels
    if (hmsh_old.nel_per_level(ilev) > 0)
      if (ncomp == 0); ncomp = size (eps_old{ilev}, 3); end
      ml = hmsh_old.msh_lev{ilev};
      ntot_old = ntot_old + ml.nel * ml.nqn;
    end
  end
  if (ncomp == 0)
    error ('cpt_transfer_hier: no active elements in the old mesh');
  end

  rdim = hmsh_old.rdim;
  Xold = zeros (ntot_old, rdim);
  Vold = zeros (ntot_old, ncomp);
  off = 0;
  for ilev = 1:hmsh_old.nlevels
    if (hmsh_old.nel_per_level(ilev) == 0); continue; end
    ml = hmsh_old.msh_lev{ilev};
    n = ml.nel * ml.nqn;
    for idim = 1:rdim
      Xold(off+(1:n), idim) = reshape (ml.geo_map(idim, :, :), n, 1);
    end
    % eps_old{ilev} is [nel x nqn x ncomp]; row order (el-major, qn-minor)
    Vold(off+(1:n), :) = reshape (permute (eps_old{ilev}, [2 1 3]), n, ncomp);
    off = off + n;
  end
  sqXold = sum (Xold.^2, 2).';                      % 1 x ntot_old

  % ---- assign each new quadrature point the nearest old value ----
  eps_new = cell (hmsh_new.nlevels, 1);
  blk = 2000;                                       % block size for the NN search
  for ilev = 1:hmsh_new.nlevels
    if (hmsh_new.nel_per_level(ilev) == 0); eps_new{ilev} = []; continue; end
    ml = hmsh_new.msh_lev{ilev};
    n = ml.nel * ml.nqn;
    Xnew = zeros (n, rdim);
    for idim = 1:rdim
      Xnew(:, idim) = reshape (ml.geo_map(idim, :, :), n, 1);
    end
    Vnew = zeros (n, ncomp);
    for a = 1:blk:n
      b = min (a + blk - 1, n);
      Xb = Xnew(a:b, :);
      % squared distances to all old points, no toolbox dependency
      D = sum (Xb.^2, 2) - 2 * (Xb * Xold.') + sqXold;
      [~, idx] = min (D, [], 2);
      Vnew(a:b, :) = Vold(idx, :);
    end
    % back to [nel x nqn x ncomp]
    eps_new{ilev} = permute (reshape (Vnew, ml.nqn, ml.nel, ncomp), [2 1 3]);
  end
end
