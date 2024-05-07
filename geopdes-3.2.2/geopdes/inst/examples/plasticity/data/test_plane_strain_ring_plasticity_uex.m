function [u_r] = test_plane_strain_ring_plasticity_uex (E, nu, sigma_y, P)

% Inner and outer radius of our ring
  a = 100;                          % inner radius mm
  b = 200;                          % outer radius mm
  Y = 2*sigma_y/sqrt(3);            % Von Mises Criterion MPa
  P_0 = Y/2*(1-(a^2)/(b^2))         % yielding pressure MPa         

  % [theta, ~] = cart2pol (x, y);
  % theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;

  if P < P_0
      P
      u_r = 2*P*b/(E*(b^2/a^2-1))*(1-nu^2);
  else
      func = @(c) log(c/a)+.5*(1-(c^2)/(b^2))-(P/Y);
      C = fsolve(func,150)
      P
      u_r = (Y*C^2)/(E*b)*(1-nu^2);
  end

end
