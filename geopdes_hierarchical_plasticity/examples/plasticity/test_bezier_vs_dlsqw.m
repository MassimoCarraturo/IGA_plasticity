% TEST_BEZIER_VS_DLSQW  Relative gap between Bezier and DLSQ_W coefficients
%
% With nquad = p+1 the element fit is square, so the fit weight cancels and the
% two operators coincide at every degree. Extra quadrature points open a gap.

function test_bezier_vs_dlsqw (varargin)

  ip = inputParser;
  ip.addParameter ('degrees', [2 3 4], @isnumeric);
  ip.addParameter ('extraquad', [0 1 2], @isnumeric);  % nquad = p+1+extra
  ip.addParameter ('nbands', 2, @isnumeric);
  ip.parse (varargin{:});
  opt = ip.Results;

  tb = projector_testbed ();  c1 = tb.setup (); %#ok<NASGU>

  f = @(x,y) 200 - abs (sqrt(x.^2+y.^2) - 160) + 20*sin (3*atan2(y,x));

  fprintf ('\n=================================================================\n');
  fprintf (' relative difference between Bezier and DLSQ_w coefficients\n');
  fprintf (' (they share the averaging; only the element fit measure differs)\n');
  fprintf ('=================================================================\n');
  fprintf ('%-8s', 'degree');
  for e = opt.extraquad
      if e == 0
          fprintf ('  %18s', 'nquad = p+1 (square)');
      else
          fprintf ('  %18s', sprintf('nquad = p+%d', 1+e));
      end
  end
  fprintf ('\n%s\n', repmat('-', 1, 8 + 20*numel(opt.extraquad)));

  for p = opt.degrees
      fprintf ('p = %-4d', p);
      for e = opt.extraquad
          d = compare (tb, p, p + 1 + e, opt.nbands, f);
          fprintf ('  %18.3e', d);
      end
      fprintf ('\n');
  end
  fprintf ('\n');
end


function d = compare (tb, p, nq, nbands, f)

  [hmsh, hspace] = tb.build (p, 5, nq, 'BEZIER', nbands, [160 40], false);
  data = tb.sample (hmsh, f);

  cB = bezier_projection_levelwise (hspace, hmsh, data);            % physical
  cD = bezier_projection_levelwise (hspace, hmsh, data, 'unit');    % identity

  d = norm (cB(:,1) - cD(:,1)) / max (norm (cB(:,1)), realmin);
end
