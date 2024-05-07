function [hspace] = MakeHspace(nurbs, truncated)
% Creation of the hspace data struct with the NURBS toolbox data structure. 
% Hspace contains information on the refinement and hierarchical basis. 
% If you are not using the NURBS Toolbox, use function CreateGeom.m instead.
%
% Input:
%   NURBS     - Data struct created in the NURBS Toolbox
% Optional input:
%   truncated - 1 for truncated refinement, 0 for normal refinement
%               Is 0 as default.
%
% Output:
%   hspace  - Hierarchical basis data structure. The cells of each value
%             correspond to the refinement level.
%
%   hspace.active{L}      - active functions of level L
%   hspace.passive{L}     - passive functions of level L
%   hspace.deactivated{L} - deactivated functions of level L
%   hspace.level          - highest level of refinement in hierarchical
%                           basis
%   hspace.dim            - dimension of the parametric space
%   hspace.truncated      - 1 if truncated and 0 for normal refinement
%   hspace.R{L1}{L2}      - Basis level operator describing basis L1 ans a
%                           linear combination of higher basis L2
%
%
% Created by Pieter van Zuijlen

 if ~exist('truncated','var')
      truncated = 0;
 end

if iscell(nurbs.knots) == 0
    dim    = 1;
    univar = 1;
    n      = length(nurbs.knots)-nurbs.order(1);
else
    dim    = numel(nurbs.knots);
    univar = 0;
    n      = 1;

    for d = 1:dim
        n = n*(length(nurbs.knots{d})-nurbs.order(d));
    end
end

assert(~any(nurbs.order-nurbs.order(1)),'The order should be the same for each dimension')
assert(dim <= 3, 'The maximum dimension of this toolbox is 3')

hspace.nfunc                 = n;
hspace.level                 = 1;
hspace.truncated             = truncated;
hspace.ref{1}.active{1}      = ones(1,n);
hspace.ref{1}.deactivated{1} = zeros(1,n);
hspace.ref{1}.passive{1}     = zeros(1,n);
hspace.R{1}{1}               = speye(n);

end