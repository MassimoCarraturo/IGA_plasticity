function M = front_kink_metrics (r, st_h, dst_h, st_e, dst_e, vm_h, sigma_y, c, dband, djump)
%FRONT_KINK_METRICS  Front-focused error metrics for a projected stress field.
%
%   M = front_kink_metrics (r, st_h, dst_h, st_e, dst_e, vm_h, sigma_y, c,
%                           dband, djump)
%
% All inputs are column vectors sampled along a ray that crosses the
% elastic-plastic front, with r the (monotone) physical distance along it:
%   st_h, dst_h   discrete hoop stress and its derivative along the ray
%   st_e, dst_e   reference hoop stress and its derivative (analytic, or an
%                 overkill numerical solution)
%   vm_h          discrete von Mises stress along the ray
%   sigma_y       yield stress (for the yield-cap violation)
%   c             front position along the ray
%   dband         half-width of the front band |r-c| <= dband
%   djump         offset at which the one-sided slopes are read, for J
%
% The rationale: the domain-wide L2 stress error integrates over a mostly
% smooth field and is nearly blind to the KINK the front leaves in the hoop
% stress (which is continuous there, while its derivative jumps).  The metrics
% below probe that feature directly.
%
% FIELDS OF M
%   E_glob   relative L2 error of st over the whole ray        (the usual metric)
%   E_band   relative L2 error of st over the front band       (localized)
%   E_inf    max |st - st_e| in the band, in % of sigma_y      (worst point)
%   E_slope  relative L2 error of d st/dr in the band          (the kink itself)
%   J        computed / reference slope jump across the front
%            (1 sharp, <1 smeared, >1 over-steepened)
%   V_yield  ||(vm - sigma_y)_+||_L2 / sigma_y over the ray    (reference-free:
%            perfect plasticity forbids vm > sigma_y, so this needs no
%            reference solution at all)
%   TV       excess total variation of st in the band, TV_h/TV_e - 1
%            (<0 the stress peak is flattened, >0 spurious oscillation)

  r = r(:);  st_h = st_h(:);  dst_h = dst_h(:);
  st_e = st_e(:);  dst_e = dst_e(:);  vm_h = vm_h(:);

  w = ones(size(r));                 % ray measure; see note below
  inband = abs(r - c) <= dband;
  all_   = true(size(r));

  rel = @(num, den, msk) sqrt(trapz(r(msk), w(msk).*num(msk).^2)) / ...
                         sqrt(trapz(r(msk), w(msk).*den(msk).^2));

  M.E_glob  = rel(st_h - st_e, st_e, all_);
  M.E_band  = rel(st_h - st_e, st_e, inband);
  M.E_inf   = 100 * max(abs(st_h(inband) - st_e(inband))) / sigma_y;
  M.E_slope = rel(dst_h - dst_e, dst_e, inband);

  % ---- slope jump across the front -------------------------------------
  rm = c - djump;  rp = c + djump;
  if rm >= min(r) && rp <= max(r)
      jh = interp1(r, dst_h, rm) - interp1(r, dst_h, rp);
      je = interp1(r, dst_e, rm) - interp1(r, dst_e, rp);
      M.J = jh / je;
  else
      M.J = NaN;
  end

  % ---- yield-cap violation (reference-free) ----------------------------
  over = max(vm_h - sigma_y, 0);
  M.V_yield = sqrt(trapz(r, over.^2)) / (sigma_y * sqrt(trapz(r, ones(size(r)))));

  % ---- excess total variation in the band ------------------------------
  tv = @(f) sum(abs(diff(f)));
  M.TV = tv(st_h(inband)) / tv(st_e(inband)) - 1;
end
