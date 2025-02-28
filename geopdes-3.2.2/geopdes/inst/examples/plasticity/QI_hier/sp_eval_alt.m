% SP_EVAL: Compute the value of a hierarchical spline function, given by its degrees of freedom, at a given set of points.%
%
% INPUT    
% u:         spline coefficients
% hspace:    hierarchical space
% uu:   (scattered) points in the parametric domain
%
% OUTPUT
% eu: values of the spline function at the given points 
% 
% Author: Cesare Bracco

function [eu] = sp_eval_alt (u, hspace, uu)

  C = hspace_subdivision_matrix (hspace, [], 'full');
  u_lev =  C{hspace.nlevels} * u;
  pl=hspace.space_of_level(1).degree;
  u_knotl=hspace.space_of_level(hspace.nlevels).knots;
  finest_basis=basisfun_multi(pl, uu, u_knotl);
  eu = u_lev'*finest_basis;

end