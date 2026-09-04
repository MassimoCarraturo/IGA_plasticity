function eps_new = wplsq_transfer_hier (hmsh_old, hspace_old, eps_old, hmsh_new)
% WPLSQ_TRANSFER_HIER  Weighted patch-based least-squares (SPR) transfer of
% quadrature-point history data between two hierarchical meshes.
%
%   eps_new = wplsq_transfer_hier (hmsh_old, hspace_old, eps_old, hmsh_new)
%
% INPUT
%   hmsh_old   : hierarchical_mesh the data currently lives on
%   hspace_old : scalar hierarchical_space on hmsh_old (its basis provides the
%                partition-of-unity blend weights)
%   eps_old    : cell over old levels, eps_old{ilev} = [nel x nqn x ncomp]
%   hmsh_new   : hierarchical_mesh to transfer the data onto
%
% OUTPUT
%   eps_new    : cell over new levels, eps_new{ilev} = [nel x nqn x ncomp]
%
% METHOD
%   Weighted patch-based least-squares fit (WPLSQ) of Hennig, Ambati,
%   De Lorenzis & Kaestner (CMAME 334, 2018, Sec. 3.3.1(iii)), the SPR
%   technique of Kumar et al. adapted to THB-splines.  Three steps:
%     1. For each old basis function N_J, the patch P_J is the set of old
%        active elements in its support, and Qbar_J the old quadrature points
%        in that patch.
%     2. A monomial field sigma_J = P a_J of the same tensor order as the
%        field basis is fitted per patch by least squares (Eq. 46), in
%        centered/scaled PHYSICAL coordinates for conditioning.  If a patch
%        holds fewer points than monomials, it is grown by adjacent elements;
%        a still-deficient patch is an error.
%     3. At each new quadrature point x_K the recovered fields are blended
%        with the old basis functions as partition-of-unity weights,
%        H[x_K] = sum_J N_J(x_K) sigma_J(x_K), which requires the truncated
%        hierarchical basis (sum_J N_J = 1).  N_J(x_K) is evaluated by
%        locating x_K in the old mesh (parametric coordinates, shared
%        geometry) and evaluating the containing element's THB shape
%        functions there.  A degree-p polynomial is reproduced exactly.

  ndim = hmsh_old.ndim;
  rdim = hmsh_old.rdim;
  if (rdim < ndim)
    error ('wplsq_transfer_hier: needs at least ndim physical coordinates');
  end
  % Monomial recovery uses the first ndim physical coordinates.  For solid
  % meshes (ring rdim=2, sphere rdim=3) this is all of them; for a planar
  % patch embedded in 3-D (rdim=3, ndim=2, constant third coordinate) it is
  % the two coordinates that vary.
  p    = hspace_old.space_of_level(1).degree;       % 1 x ndim
  nmono = prod (p + 1);                             % tensor monomials up to p

  % ==== gather old active elements: coords, values, global THB dofs ========
  ncomp = 0;
  for ilev = 1:hmsh_old.nlevels
    if (hmsh_old.nel_per_level(ilev) > 0)
      ncomp = size (eps_old{ilev}, 3); break
    end
  end
  if (ncomp == 0); error ('wplsq_transfer_hier: no active elements in the old mesh'); end

  nel_tot = sum (hmsh_old.nel_per_level);
  E = repmat (struct ('lev', [], 'Xphys', [], 'val', [], 'dofs', []), nel_tot, 1);
  dof2elem = cell (hspace_old.ndof, 1);
  ie = 0;
  for ilev = 1:hmsh_old.nlevels
    if (hmsh_old.nel_per_level(ilev) == 0); continue; end
    ml   = hmsh_old.msh_lev{ilev};
    sp   = sp_evaluate_element_list (hspace_old.space_of_level(ilev), ml, 'value', true);
    CsubT = hspace_old.Csub{ilev}.';
    for iel = 1:ml.nel
      ie = ie + 1;
      conn = sp.connectivity(:, iel); mask = conn > 0; conn = conn(mask);
      dofs = find (any (CsubT(:, conn), 2));         % global THB dofs on element
      E(ie).lev   = ilev;
      E(ie).Xphys = reshape (ml.geo_map(:, :, iel), rdim, ml.nqn).';   % nqn x rdim
      E(ie).val   = reshape (eps_old{ilev}(iel, :, :), ml.nqn, ncomp); % nqn x ncomp
      E(ie).dofs  = dofs(:).';
      for J = E(ie).dofs
        dof2elem{J}(end+1) = ie;                     %#ok<AGROW>
      end
    end
  end
  nel_tot = ie;

  % ==== per-patch monomial recovery ========================================
  % Store, per dof J that is ever needed, its fit (coeffs, center, scale).
  fitA      = cell (hspace_old.ndof, 1);
  fitCenter = cell (hspace_old.ndof, 1);
  fitScale  = cell (hspace_old.ndof, 1);

  for J = 1:hspace_old.ndof
    els = dof2elem{J};
    if (isempty (els)); continue; end
    % grow the patch until it holds enough points (Hennig boundary rule)
    npts = patch_npts (E, els);
    grow = 0;
    while (npts < nmono && grow < 3)
      els = grow_patch (E, dof2elem, els);
      npts = patch_npts (E, els);
      grow = grow + 1;
    end
    X = vertcat (E(els).Xphys);                      % npts x rdim
    F = vertcat (E(els).val);                        % npts x ncomp
    if (size (X, 1) < nmono)
      error (['wplsq_transfer_hier: patch of dof %d has %d points for %d ' ...
              'monomials even after growth'], J, size (X, 1), nmono);
    end
    c = mean (X, 1);
    s = max (X, [], 1) - min (X, [], 1);
    s(s == 0) = 1;
    P = build_mono (X, p, c, s, ndim);               % npts x nmono
    fitA{J}      = P \ F;                             % nmono x ncomp (LSQ)
    fitCenter{J} = c;
    fitScale{J}  = s;
  end

  % ==== blend at the new quadrature points =================================
  % precompute level break points and nel_dir for point location
  brk = cell (hmsh_old.nlevels, 1);
  neld = cell (hmsh_old.nlevels, 1);
  for l = 1:hmsh_old.nlevels
    brk{l}  = hmsh_old.mesh_of_level(l).breaks;
    neld{l} = hmsh_old.mesh_of_level(l).nel_dir;
  end

  eps_new = cell (hmsh_new.nlevels, 1);
  for ilev = 1:hmsh_new.nlevels
    if (hmsh_new.nel_per_level(ilev) == 0); eps_new{ilev} = []; continue; end
    ml = hmsh_new.msh_lev{ilev};
    qnp = my_msh_evaluate_qn (hmsh_new.mesh_of_level(ilev), hmsh_new.active{ilev});  % ndim x nqn x nel parametric
    out = zeros (ml.nel, ml.nqn, ncomp);
    for iel = 1:ml.nel
      for iqn = 1:ml.nqn
        xi = qnp(:, iqn, iel).';                      % 1 x ndim parametric
        xp = ml.geo_map(:, iqn, iel).';              % 1 x rdim physical
        lev = locate_level (hmsh_old, brk, neld, xi);
        [gdofs, w] = thb_at_param (hspace_old, lev, xi, ndim);
        acc = zeros (1, ncomp);
        wsum = 0;
        for k = 1:numel (gdofs)
          J = gdofs(k);
          if (w(k) == 0 || isempty (fitA{J})); continue; end
          z = build_mono (xp, p, fitCenter{J}, fitScale{J}, ndim);
          acc = acc + w(k) * (z * fitA{J});
          wsum = wsum + w(k);
        end
        if (wsum > 0); acc = acc / wsum; end         % guard partial PU at edges
        out(iel, iqn, :) = acc;
      end
    end
    eps_new{ilev} = out;
  end
end

% ------------------------------------------------------------------------
function n = patch_npts (E, els)
  n = 0;
  for e = els; n = n + size (E(e).Xphys, 1); end
end

function els = grow_patch (E, dof2elem, els)
% add every element sharing a dof with the current patch (one ring)
  dofs = unique ([E(els).dofs]);
  add = [];
  for J = dofs; add = [add, dof2elem{J}]; end        %#ok<AGROW>
  els = unique ([els, add]);
end

function P = build_mono (X, p, c, s, ndim)
% tensor monomials z1^i z2^j (z3^k), 0<=.<=p, z = (X-c)./s ; X is n x ndim
  z = (X - c) ./ s;
  pw = cell (1, ndim);
  for d = 1:ndim
    pw{d} = z(:, d) .^ (0:p(d));                     % n x (p(d)+1)
  end
  P = pw{1};
  for d = 2:ndim
    nP = size (P, 2); nD = size (pw{d}, 2);
    P = reshape (P, [], nP, 1) .* reshape (pw{d}, [], 1, nD);
    P = reshape (P, size (z, 1), nP * nD);
  end
end

function lev = locate_level (hmsh, brk, neld, xi)
% finest active level whose cell contains parametric point xi
  lev = 0;
  for l = hmsh.nlevels:-1:1
    sub = zeros (1, numel (xi));
    ok = true;
    for d = 1:numel (xi)
      b = brk{l}{d};
      k = find (xi(d) >= b(1:end-1) - 1e-12 & xi(d) <= b(2:end) + 1e-12, 1, 'last');
      if (isempty (k)); ok = false; break; end
      sub(d) = k;
    end
    if (~ok); continue; end
    if (numel (xi) == 1)
      lin = sub(1);
    else
      subc = num2cell (sub);
      lin = sub2ind (neld{l}, subc{:});
    end
    if (any (hmsh.active{l} == lin)); lev = l; return; end
  end
  if (lev == 0)
    error ('wplsq_transfer_hier: could not locate a new quadrature point in the old mesh');
  end
end

function [gdofs, vals] = thb_at_param (hspace, lev, xi, ndim)
% THB shape functions of hspace at parametric point xi (in an active level-lev
% element): level-lev TP B-splines at xi, mapped through Csub{lev}
  sp = hspace.space_of_level(lev);
  gidx = cell (1, ndim); gval = cell (1, ndim);
  for d = 1:ndim
    U = sp.knots{d}; pd = sp.degree(d); nd = sp.ndof_dir(d);
    span = findspan (nd - 1, pd, xi(d), U);
    B = basisfun (span, xi(d), pd, U);               % 1 x (pd+1)
    gidx{d} = (span - pd + 1) : (span + 1);          % 1-based global TP indices
    gval{d} = B(:).';
  end
  subs = cell (1, ndim); vsub = cell (1, ndim);
  [subs{1:ndim}] = ndgrid (gidx{:});
  [vsub{1:ndim}] = ndgrid (gval{:});
  lin = sub2ind (sp.ndof_dir, subs{:});
  vv = vsub{1};
  for d = 2:ndim; vv = vv .* vsub{d}; end
  lin = lin(:); vv = vv(:);
  thb = (vv.' * hspace.Csub{lev}(lin, :)).';         % last_dof(lev) x 1
  gdofs = find (abs (thb) > 1e-13).';
  vals  = thb(gdofs).';
end
