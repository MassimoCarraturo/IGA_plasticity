% TEST_PROJECTOR_CONVERGENCE  Relative L2 error of a single projection of a known function under mesh refinement

function test_projector_convergence (varargin)

  ip = inputParser;
  ip.addParameter ('fun', 'smooth', @ischar);
  ip.addParameter ('family', 'uniform', @ischar);
  ip.addParameter ('p', 4, @isnumeric);
  ip.addParameter ('proj', {'L2','QI','QI_C0','Bezier','Bezier_C0','DLSQ','DLSQ_W'}, @iscell);
  ip.addParameter ('lambda', [], @isnumeric);   % QI Tikhonov override
  ip.parse (varargin{:});
  opt = ip.Results;

  tb = projector_testbed ();  c1 = tb.setup (); %#ok<NASGU>

  f = test_function (opt.fun);
  switch lower (opt.family)
      case 'uniform';      steps = {5, 10, 20, 40};
      case 'hierarchical'; steps = {0, 1, 2, 3};
      otherwise; error ('family must be uniform or hierarchical');
  end

  fprintf ('\n==================================================================\n');
  fprintf (' projector convergence, function = %s, family = %s, p = %d\n', ...
           opt.fun, opt.family, opt.p);
  fprintf ('==================================================================\n');

  nS = numel (steps);
  err = nan (numel (opt.proj), nS);
  ndof = nan (1, nS);  nelv = nan (1, nS);

  for is = 1:nS
      for ipj = 1:numel (opt.proj)
          [e, nd, ne] = one_case (tb, opt.p, opt.family, steps{is}, opt.proj{ipj}, f, opt.lambda);
          err(ipj, is) = e;  ndof(is) = nd;  nelv(is) = ne;
      end
      fprintf ('  mesh %d: %d elements, %d dof\n', is, nelv(is), ndof(is));
  end

  fprintf ('\n%-12s', 'operator');
  for is = 1:nS; fprintf ('  %11s', sprintf('nel=%d', nelv(is))); end
  fprintf ('     rate    verdict\n');
  fprintf ('%s\n', repmat('-', 1, 78));
  for ipj = 1:numel (opt.proj)
      fprintf ('%-12s', opt.proj{ipj});
      for is = 1:nS; fprintf ('  %11.3e', err(ipj, is)); end
      e = err(ipj, :);
      ok = all (isfinite (e));
      if ok && all (e < 1e-11)
          fprintf ('        --   EXACT\n');
      elseif ok
          r = log (e(1)/e(end)) / log (nelv(end)/nelv(1)) * 2;   % rate in h
          mono = all (e(2:end) <= e(1:end-1) * (1 + 1e-12));
          if mono; v = 'converges'; else; v = '*** NOT MONOTONE ***'; end
          fprintf ('  %6.2f   %s\n', r, v);
      else
          fprintf ('        --   FAILED\n');
      end
  end
  fprintf ('\n');
end


function [err, ndof, nel] = one_case (tb, p, family, step, proj, f, lam)

  if strcmpi (family, 'uniform')
      nsub = step;  nbands = 0;
  else
      nsub = 5;     nbands = step;
  end
  % Band width must stay fixed, a band that halves each step never contains a full (p+1)^d support and no function activates
  [hmsh, hspace] = tb.build (p, nsub, p+1, proj, nbands, [160 40], true);

  eps_pl = tb.sample (hmsh, f);
  coef = history_variable_projection_hier (hspace, hmsh, eps_pl, proj, hmsh, lam);
  back = tb.evaluate (hmsh, hspace, coef);

  num = 0;  den = 0;
  for ilev = 1:hmsh.nlevels
      if hmsh.nel_per_level(ilev) == 0; continue; end
      ml = hmsh.msh_lev{ilev};
      w  = ml.quad_weights .* ml.jacdet;
      for iel = 1:ml.nel
          val = back{ilev}(iel,:,1).';  ex = eps_pl{ilev}(iel,:,1).';
          num = num + sum (w(:,iel) .* (val - ex).^2);
          den = den + sum (w(:,iel) .* ex.^2);
      end
  end
  err = sqrt (num / den);
  ndof = hspace.ndof;  nel = hmsh.nel;
end


function f = test_function (name)
  switch lower (name)
      case 'poly'                       % degree 4, must be reproduced exactly
          f = @(x,y) 1 + 0.01*x - 0.02*y + 1e-4*x.*y + 1e-6*x.^2.*y.^2;
      case 'smooth'                     % analytic, rate p+1 expected
          f = @(x,y) 200 * sin (pi*sqrt(x.^2+y.^2)/200) .* cos (2*atan2(y,x));
      case 'kink'                       % C^0 ridge at r = 160, mimics the front
          f = @(x,y) 200 - abs (sqrt(x.^2+y.^2) - 160);
      otherwise
          error ('unknown test function %s', name);
  end
end
