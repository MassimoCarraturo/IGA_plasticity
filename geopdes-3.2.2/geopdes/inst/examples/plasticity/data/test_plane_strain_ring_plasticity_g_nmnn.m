function g = test_plane_strain_ring_plasticity_g_nmnn (x, y, P, nu, ~)

% Inner and outer radius of our ring
  R_i = 100; R_o = 200;

  g = zeros(2, size(x,1), size(x,2));
  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;

  %C = (P * R_i^2) / (R_o^2 - R_i^2);
  %C = C*((1-nu)/((1+nu)*(1-2*nu)) - 1);
    C = P;
   g(1, :, :) = C * cos (theta);
   g(2, :, :) = C * sin (theta);
end

