% TEST_PROJECTOR_REPRODUCTION  Check that each projector reproduces a spline lying in its own hierarchical space
% Speleers and Manni, Numer. Math. 132 (2016): the level-wise QI against the truncated basis reproduces the hierarchical space

function test_projector_reproduction (varargin)

  ip = inputParser;
  ip.addParameter ('p', 4, @isnumeric);
  ip.addParameter ('nsub', 5, @isnumeric);
  ip.addParameter ('nbands', [0 1 2], @isnumeric);
  ip.addParameter ('proj', {'L2','QI','QI_C0','Bezier','Bezier_C0','DLSQ','DLSQ_W'}, @iscell);
  ip.addParameter ('lambda', 0, @isnumeric);
  ip.parse (varargin{:});
  opt = ip.Results;

  tb = projector_testbed ();  c1 = tb.setup (); %#ok<NASGU>

  fprintf ('\n===============================================================\n');
  fprintf (' reproduction of the operator''s own space, p = %d\n', opt.p);
  fprintf (' relative max error of the returned coefficients\n');
  fprintf ('===============================================================\n');
  fprintf ('%-12s', 'operator');
  for nb = opt.nbands; fprintf ('  %13s', sprintf('%d level(s)', nb+1)); end
  fprintf ('\n%s\n', repmat('-', 1, 12 + 15*numel(opt.nbands)));

  for ipj = 1:numel (opt.proj)
      fprintf ('%-12s', opt.proj{ipj});
      for nb = opt.nbands
          e = one_case (tb, opt, opt.proj{ipj}, nb);
          if (isnan (e))
              fprintf ('  %13s', 'n/a');
          else
              fprintf ('  %13.3e', e);
          end
      end
      fprintf ('\n');
  end
  fprintf ('\n');
end


function err = one_case (tb, opt, proj, nbands)

  [hmsh, hspace] = tb.build (opt.p, opt.nsub, opt.p+1, proj, nbands, [160 40], false);

  rng (7);
  c0 = 1 + 0.5 * sin ((1:hspace.ndof).');

  data = tb.evaluate (hmsh, hspace, c0);
  coef = history_variable_projection_hier (hspace, hmsh, data, proj, hmsh, opt.lambda);
  c1 = coef(:, 1);

  if (~all (isfinite (c1)))
      err = nan;  return
  end
  % Field L2 error rather than coefficient max norm, a truncated function with tiny support can carry a large coefficient error while the field barely moves
  err = tb.rel_err (hmsh, tb.evaluate (hmsh, hspace, c1), data);
end
