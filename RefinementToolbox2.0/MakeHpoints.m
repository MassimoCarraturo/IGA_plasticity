function [hpoints] = MakeHpoints(nurbs,dim)
% This struct contains the control points for the hierarcical basis.
% Different cells contain the points of different refinements.
%
% Input:
%   NURBS - Data struct created in the NURBS Toolbox
%   dim   - Dimension of the physical space
%
% Output:
%   hpoints  - Cells with control points
%
%
%
% Created by Pieter van Zuijlen

cpts = reshape(nurbs.coefs,[4,prod(nurbs.number)]);

if ~exist('dim','var')
    dim = size(cpts,1);
end

hpoints{1} = cpts(1:dim,:)';

end