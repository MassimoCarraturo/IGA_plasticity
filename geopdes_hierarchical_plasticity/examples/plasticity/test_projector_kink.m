% TEST_PROJECTOR_KINK  Realised stability constant of each projector on a kinked field
% amp = err / err of the global L2 fit in the same target space, since ||f - Pf|| <= ||P|| dist(f,V)

function test_projector_kink (varargin)

  ip = inputParser;
  ip.addParameter ('p', 4, @isnumeric);
  ip.addParameter ('nsub', 5, @isnumeric);
  ip.addParameter ('nbands', [0 1 2], @isnumeric);
  ip.addParameter ('rfront', 160, @isnumeric);   % 160 sits ON an element boundary
  ip.addParameter ('bandc', 160, @isnumeric);    % centre of the refined band
  ip.addParameter ('bandw', 40, @isnumeric);     % its half-width
  ip.addParameter ('proj', {'L2','QI','Bezier','DLSQ','L2_C0','QI_C0','Bezier_C0'}, @iscell);
  ip.parse (varargin{:});
  opt = ip.Results;

  tb = projector_testbed ();  c1 = tb.setup (); %#ok<NASGU>

  rfront = opt.rfront;
  fkink = @(x,y) max (0, rfront - sqrt (x.^2 + y.^2));

  fprintf ('\n');
  fprintf ('=========================================================================\n');
  fprintf (' one projection of a kinked field, p = %d, ring, front at r = %g\n', opt.p, rfront);
  fprintf (' err  = relative L2 error of the projected field\n');
  fprintf (' amp  = err / err of the BEST approximation in the same target space\n');
  fprintf ('=========================================================================\n');

  for nb = opt.nbands
    fprintf ('\n--- %d level(s) ---\n', nb + 1);
    fprintf ('%-12s %8s %14s %10s\n', 'operator', 'ndof', 'err', 'amp');
    fprintf ('%s\n', repmat ('-', 1, 48));

    e = struct ();  n = struct ();
    for ipj = 1:numel (opt.proj)
      tag = opt.proj{ipj};
      [ee, nn] = one_case (tb, opt, tag, nb, fkink);
      e.(tag) = ee;  n.(tag) = nn;
    end

    for ipj = 1:numel (opt.proj)
      tag = opt.proj{ipj};
      base = base_of (tag);
      if (isnan (e.(tag)))
        fprintf ('%-12s %8d %14s %10s\n', tag, n.(tag), 'diverged', '-');
      else
        fprintf ('%-12s %8d %14.4e %10.2f\n', tag, n.(tag), e.(tag), e.(tag) / e.(base));
      end
    end

    if (isfield (e, 'L2') && isfield (e, 'L2_C0'))
      fprintf ('\n  space effect at fixed rule (global L2):  err(C^{p-1}) / err(C^0) = %.2f\n', ...
               e.L2 / e.L2_C0);
      fprintf ('  and the C^0 space carries %.1f times the functions\n', n.L2_C0 / n.L2);
    end
  end
  fprintf ('\n');
end


function b = base_of (tag)
  if (numel (tag) > 3 && strcmp (tag(end-2:end), '_C0'))
    b = 'L2_C0';
  else
    b = 'L2';
  end
end


function [err, ndof] = one_case (tb, opt, tag, nbands, fkink)

  % L2_C0 is the global L2 fit on the C0 space, the reference for the C0 operators
  if (strcmp (tag, 'L2_C0')); proj = 'L2'; else; proj = tag; end

  [hmsh, hspace] = tb.build (opt.p, opt.nsub, opt.p+1, tag, nbands, [opt.bandc opt.bandw], false);
  ndof = hspace.ndof;

  data = tb.sample (hmsh, fkink);
  coef = history_variable_projection_hier (hspace, hmsh, data, proj, hmsh, 0);
  if (~all (isfinite (coef(:, 1))))
    err = nan;  return
  end
  err = tb.rel_err (hmsh, tb.evaluate (hmsh, hspace, coef(:, 1)), data);
end
