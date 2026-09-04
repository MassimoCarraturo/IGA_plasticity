% TEST_PROJECTOR_STABILITY  Drift of a transfer operator under repeated project-evaluate cycles
% mode 'fixed' iterates on one mesh, 'cycle' refines and coarsens between projections

function test_projector_stability (varargin)

  ip = inputParser;
  ip.addParameter ('mode', 'fixed', @ischar);
  ip.addParameter ('niter', 12, @isnumeric);
  ip.addParameter ('p', 4, @isnumeric);
  ip.addParameter ('nsub', 5, @isnumeric);
  ip.addParameter ('nbands', 2, @isnumeric);
  ip.addParameter ('proj', {'L2','QI','QI_C0','Bezier','Bezier_C0','DLSQ','DLSQ_W'}, @iscell);
  ip.addParameter ('lambda', [], @isnumeric);
  ip.parse (varargin{:});
  opt = ip.Results;

  tb = projector_testbed ();  c1 = tb.setup (); %#ok<NASGU>

  f = @(x,y) 200 - abs (sqrt(x.^2+y.^2) - 160) + 20*sin(3*atan2(y,x));

  fprintf ('\n=================================================================\n');
  fprintf (' repeated application, mode = %s, %d iterations, p = %d\n', ...
           opt.mode, opt.niter, opt.p);
  fprintf ('=================================================================\n');
  fprintf ('%-12s %10s %10s %10s %10s   %s\n', ...
           'operator', 'err(1)', 'err(2)', sprintf('err(%d)',opt.niter), ...
           'max|c|', 'verdict');
  fprintf ('%s\n', repmat('-', 1, 76));

  for ipj = 1:numel (opt.proj)
      [e, cmax] = run_one (tb, opt, opt.proj{ipj}, f);
      if any (~isfinite (e))
          fprintf ('%-12s %10s %10s %10s %10s   DIVERGED\n', opt.proj{ipj}, '-','-','-','-');
          continue
      end
      growth = e(end) / max (e(2), realmin);
      if growth > 10
          v = sprintf ('*** GROWS x%.3g ***', growth);
      elseif growth > 1.05
          v = sprintf ('drifts x%.2f', growth);
      else
          v = 'stable';
      end
      fprintf ('%-12s %10.3e %10.3e %10.3e %10.3e   %s\n', ...
               opt.proj{ipj}, e(1), e(2), e(end), cmax, v);
  end
  fprintf ('\n');
end


function [err, cmax] = run_one (tb, opt, proj, f)

  [hmshA, hspA] = tb.build (opt.p, opt.nsub, opt.p+1, proj, opt.nbands, [160 40], false);
  ad = tb.refine_opts (opt.p, opt.nbands + 1, true);

  data = tb.sample (hmshA, f);
  ref  = data;                       % the exact field, for the error norm
  err  = nan (1, opt.niter);
  cmax = 0;

  for it = 1:opt.niter
      coef = history_variable_projection_hier (hspA, hmshA, data, proj, hmshA, opt.lambda);
      cmax = max (cmax, max (abs (coef(:,1))));
      if (~all (isfinite (coef(:))))
          err(it:end) = inf;  return
      end

      if strcmpi (opt.mode, 'cycle')
          % Cref and Ccoar are exact between nested spaces, so any drift is the projector and the quadrature round trip
          mk = tb.mark_band (hmshA, 160, 40);
          [hmshB, hspB, Cref] = adaptivity_refine (hmshA, hspA, mk, ad);
          dB = tb.evaluate (hmshB, hspB, Cref * coef);
          cB = history_variable_projection_hier (hspB, hmshB, dB, proj, hmshB, opt.lambda);
          cmax = max (cmax, max (abs (cB(:,1))));

          mkc = mark_finest (hmshB);
          [hmshA, hspA, Ccoar] = adaptivity_coarsen (hmshB, hspB, mkc, ad);
          data = tb.evaluate (hmshA, hspA, Ccoar * cB);
          ref  = tb.sample (hmshA, f);     % the mesh may not return identical
      else
          data = tb.evaluate (hmshA, hspA, coef);
      end

      err(it) = tb.rel_err (hmshA, data, ref);
  end
end


function marked = mark_finest (hmsh)
% Finest non-empty level, so coarsening undoes the last refinement
  marked = cell (hmsh.nlevels, 1);
  for ilev = 1:hmsh.nlevels; marked{ilev} = []; end
  for ilev = hmsh.nlevels:-1:1
      if hmsh.nel_per_level(ilev) > 0
          marked{ilev} = hmsh.active{ilev}(:).';
          return
      end
  end
end
