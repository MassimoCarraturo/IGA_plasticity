% Verification of BEZIER_PROJECTION_LEVELWISE against the properties a
% projection operator must have, on a hierarchical mesh, at several degrees.
%
% Test 1 (reproduction).  On an affine geometry a polynomial of degree <= p lies
% in the spline space, so ANY projection must return it exactly.  This is the
% test that fails for a rank-deficient element-local fit, and it is checked at
% p = 2, 3 and 4 with two and three levels.
%
% Test 2 (agreement at p = 2).  Where the element-local operator is well posed
% the two should agree closely; a large discrepancy would mean the level-wise
% variant is a different operator rather than an extension.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));
project_root = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(project_root,'nurbs-1.4.3','nurbs-1.4.3','inst')));
addpath(genpath(fullfile(project_root,'geopdes-3.2.2','geopdes','inst')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','hierarchical_classes')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','adaptivity_iga')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','initialize')));
addpath(genpath(fullfile(project_root,'geopdes_hierarchical_plasticity','examples','plasticity')));
cd(here);
warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix');

fprintf('%-4s %-7s %-7s %-8s %-14s %-14s\n', ...
        'p','levels','ndof','nel','LW reprod.','elem-local');

for p = [2 3 4]
for nref = [1 2]

  % ---- unit square, identity map: physical = parametric, so polynomials of
  % degree <= p are exactly representable
  geometry = geo_load ('geo_square.txt');
  [knots, zeta] = kntrefine (geometry.nurbs.knots, [3 3], [p p], [p-1 p-1]);
  rule     = msh_gauss_nodes ([p+1 p+1]);
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
  msh      = msh_cartesian (zeta, qn, qw, geometry);
  space    = sp_bspline (knots, [p p], msh);

  hmsh   = hierarchical_mesh (msh, [2 2]);
  hspace = hierarchical_space (hmsh, space, 'standard', true, [p-1 p-1]);

  % ---- refine a corner region nref times to build the hierarchy
  for k = 1:nref
    xg = zeros (hmsh.nel, 1); yg = zeros (hmsh.nel, 1);
    ct = 0;
    for ilev = 1:hmsh.nlevels
      if (hmsh.nel_per_level(ilev) == 0); continue; end
      ml = hmsh.msh_lev{ilev};
      cx = squeeze (mean (ml.geo_map(1,:,:), 2));
      cy = squeeze (mean (ml.geo_map(2,:,:), 2));
      xg(ct+(1:ml.nel)) = cx(:); yg(ct+(1:ml.nel)) = cy(:);
      ct = ct + ml.nel;
    end
    marked_all = find (xg < 0.5 & yg < 0.5);
    marked = cell (hmsh.nlevels, 1);
    off = 0;
    for ilev = 1:hmsh.nlevels
      n = hmsh.nel_per_level(ilev);
      sel = marked_all(marked_all > off & marked_all <= off + n) - off;
      marked{ilev} = hmsh.active{ilev}(sel);
      off = off + n;
    end
    adaptivity_data.flag = 'elements';
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
  end

  % ---- sample a polynomial of degree exactly p at the quadrature points
  fpoly = @(x,y) 1 + 2*x + 3*y - x.*y + x.^p + 0.5*y.^p;
  data = cell (hmsh.nlevels, 1);
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) == 0)
      data{ilev} = [];
      continue
    end
    ml = hmsh.msh_lev{ilev};
    X = permute (ml.geo_map(1,:,:), [3 2 1]);          % nel x nqn
    Y = permute (ml.geo_map(2,:,:), [3 2 1]);
    data{ilev} = reshape (fpoly (X, Y), ml.nel, ml.nqn, 1);
  end

  % ---- level-wise operator
  err_lw = NaN; err_el = NaN;
  try
    c_lw = bezier_projection_levelwise (hspace, hmsh, data);
    err_lw = reproduction_error (c_lw, hspace, hmsh, fpoly);
  catch ME
    err_lw = -1; fprintf('   LW  ERR: %s\n', ME.message);
  end
  try
    c_el = bezier_projection_hier (hspace, hmsh, data);
    err_el = reproduction_error (c_el, hspace, hmsh, fpoly);
  catch ME
    err_el = -1; fprintf('   EL  ERR: %s\n', ME.message);
  end

  s_lw = fmt (err_lw); s_el = fmt (err_el);
  fprintf('%-4d %-7d %-7d %-8d %-14s %-14s\n', p, hmsh.nlevels, hspace.ndof, hmsh.nel, s_lw, s_el);
end
end

function s = fmt (e)
  if (e < 0); s = 'FAILED'; else; s = sprintf('%.3e', e); end
end

function err = reproduction_error (coeff, hspace, hmsh, fpoly)
% relative L2 error between the projected spline and the exact polynomial
  n2 = 0; d2 = 0;
  for ilev = 1:hmsh.nlevels
    if (hmsh.nel_per_level(ilev) == 0); continue; end
    ml = hmsh.msh_lev{ilev};
    sp = sp_evaluate_element_list (hspace.space_of_level(ilev), ml, 'value', true);
    Cs = hspace.Csub{ilev};
    for iel = 1:ml.nel
      conn = sp.connectivity(:, iel); mask = conn > 0; conn = conn(mask);
      Ntp  = sp.shape_functions(:, mask, iel);
      E    = full (Cs(conn, :));                        % level-l TP x active-so-far THB
      nf   = size (E, 2);                               % Csub has last_dof(l) columns
      vals = Ntp * (E * coeff(1:nf, :));
      X = ml.geo_map(1,:,iel).'; Y = ml.geo_map(2,:,iel).';
      ex = fpoly (X, Y);
      wq = ml.quad_weights(:, iel) .* ml.jacdet(:, iel);
      n2 = n2 + sum (wq .* (vals - ex).^2);
      d2 = d2 + sum (wq .* ex.^2);
    end
  end
  err = sqrt (n2 / d2);
end
