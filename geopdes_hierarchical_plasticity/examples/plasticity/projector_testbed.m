% PROJECTOR_TESTBED  Pieces shared by the test_projector_* scripts
% tb = projector_testbed () returns handles to setup, build, refine_opts, mark_band, sample, evaluate, rel_err

function tb = projector_testbed ()
  tb.setup = @setup;  tb.build = @build;  tb.refine_opts = @refine_opts;  tb.mark_band = @mark_band;
  tb.sample = @sample;  tb.evaluate = @evaluate;  tb.rel_err = @rel_err;
end


function c = setup ()
% Paths, warnings and the working directory, restored when the returned onCleanup dies
  here = fileparts (mfilename ('fullpath'));
  root = fullfile (here, '..', '..', '..');
  addpath (genpath (fullfile (root, 'nurbs-1.4.3', 'nurbs-1.4.3', 'inst')));
  addpath (genpath (fullfile (root, 'geopdes-3.2.2', 'geopdes', 'inst')));
  addpath (genpath (fullfile (root, 'geopdes_hierarchical_plasticity', 'hierarchical_classes')));
  addpath (genpath (fullfile (root, 'geopdes_hierarchical_plasticity', 'adaptivity_iga')));
  addpath (genpath (fullfile (root, 'geopdes_hierarchical_plasticity', 'initialize')));
  addpath (genpath (fullfile (root, 'geopdes_hierarchical_plasticity', 'examples', 'plasticity')));
  addpath (genpath (fullfile (root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical')));
  addpath (fullfile (root, 'geopdes_hierarchical_plasticity', 'quasi_interpolant_hierarchical', 'libqi', 'matlab'));

  old = cd (here);  c = onCleanup (@() cd (old));
  warning ('off', 'MATLAB:nearlySingularMatrix');
  warning ('off', 'MATLAB:singularMatrix');
end


function [hmsh, hspace] = build (p, nsub, nquad, proj, nbands, band, strict)
% Ring, nsub^2 coarse elements, nbands dyadic refinements of the band |r - band(1)| <= band(2)
% The *_C0 operators get the C0 space, strict errors when a new level carries no active function
  pd.geo_name = 'geo_ring_SouzaNeto.txt';
  md.degree = [p p];  md.regularity = [p-1 p-1];
  md.nsub_coarse = [nsub nsub];  md.nsub_refine = [2 2];
  md.nquad = [nquad nquad];  md.space_type = 'standard';  md.truncated = 1;

  [hmsh, ~, geometry] = adaptivity_initialize_vector (pd, md);

  if any (strcmpi (proj, {'QI_C0', 'BEZIER_C0', 'L2_C0'}))
      reg = zeros (1, 2);
  else
      reg = md.regularity;
  end
  [knots, zeta] = kntrefine (geometry.nurbs.knots, md.nsub_coarse-1, md.degree, reg);
  rule = msh_gauss_nodes (md.nquad);
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
  msh = msh_cartesian (zeta, qn, qw, geometry);
  sp  = sp_bspline (knots, md.degree, msh);
  hspace = hierarchical_space (hmsh, sp, md.space_type, md.truncated, reg);

  ad = refine_opts (p, 1 + nbands, false);
  for k = 1:nbands
      marked = mark_band (hmsh, band(1), band(2));
      [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, ad);
      if strict && hspace.ndof_per_level(end) == 0
          error (['projector_testbed: the new level carries no active ' ...
                  'functions, so the space cannot gain resolution. Widen the band.']);
      end
  end
end


function ad = refine_opts (p, max_level, coarsen)
% adaptivity_refine needs the full adaptivity_data struct, a partial one refines the mesh but not the space
  ad.flag = 'elements';
  ad.adm_strategy = 'admissible';  ad.adm = p;
  ad.max_level = max_level;
  ad.mark_strategy = 'MS';  ad.mark_param = 0.5;
  if coarsen
      ad.mark_param_coarsening = 0.5;  ad.coarsening_flag = 'any';
  else
      ad.mark_param_coarsening = 0;    ad.coarsening_flag = 'none';
  end
  ad.max_ndof = 1e7;  ad.max_nel = 1e6;  ad.num_max_iter = 1;  ad.tol = 1e-10;
end


function marked = mark_band (hmsh, rc, halfwidth)
% Active elements whose centre lies within halfwidth of radius rc
  marked = cell (hmsh.nlevels, 1);
  for ilev = 1:hmsh.nlevels
      marked{ilev} = [];
      if hmsh.nel_per_level(ilev) == 0; continue; end
      ml = hmsh.msh_lev{ilev};
      keep = false (ml.nel, 1);
      for iel = 1:ml.nel
          g = ml.geo_map(:,:,iel);
          r = mean (sqrt (g(1,:).^2 + g(2,:).^2));
          keep(iel) = abs (r - rc) <= halfwidth;
      end
      marked{ilev} = hmsh.active{ilev}(keep);
  end
end


function data = sample (hmsh, f)
% f(x,y) at the quadrature points, in the first component of a history-variable cell
  data = cell (hmsh.nlevels, 1);
  for ilev = 1:hmsh.nlevels
      if hmsh.nel_per_level(ilev) == 0
          data{ilev} = zeros (0, hmsh.mesh_of_level(ilev).nqn, 6);  continue
      end
      ml = hmsh.msh_lev{ilev};
      v = zeros (ml.nel, ml.nqn, 6);
      for iel = 1:ml.nel
          g = ml.geo_map(:,:,iel);
          v(iel,:,1) = f (g(1,:), g(2,:));
      end
      data{ilev} = v;
  end
end


function data = evaluate (hmsh, hspace, coef)
% The hierarchical function with coefficients coef(:,1) at the quadrature points, same layout as sample
  data = cell (hmsh.nlevels, 1);
  for ilev = 1:hmsh.nlevels
      if hmsh.nel_per_level(ilev) == 0
          data{ilev} = zeros (0, hmsh.mesh_of_level(ilev).nqn, 6);  continue
      end
      ml  = hmsh.msh_lev{ilev};
      spl = sp_evaluate_element_list (hspace.space_of_level(ilev), ml, 'value', true);
      Cs  = hspace.Csub{ilev};
      ctp = Cs * coef(1:size(Cs,2), 1);
      v = zeros (ml.nel, ml.nqn, 6);
      for iel = 1:ml.nel
          cn = spl.connectivity(:, iel);  m = cn > 0;
          v(iel,:,1) = spl.shape_functions(:, m, iel) * ctp(cn(m));
      end
      data{ilev} = v;
  end
end


function e = rel_err (hmsh, a, b)
% Relative L2 distance between two sampled fields, b is the reference
  num = 0;  den = 0;
  for ilev = 1:hmsh.nlevels
      if hmsh.nel_per_level(ilev) == 0; continue; end
      ml = hmsh.msh_lev{ilev};
      w  = ml.quad_weights .* ml.jacdet;
      da = a{ilev}(:,:,1).';  db = b{ilev}(:,:,1).';
      num = num + sum (sum (w .* (da - db).^2));
      den = den + sum (sum (w .* db.^2));
  end
  e = sqrt (num / den);
end
