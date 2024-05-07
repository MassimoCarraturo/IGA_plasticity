function [X,Y,Z,hfuncs,shape] = ConstructGeometry(basis,hspace,hpoints)
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

if ~exist('nref','var')
    nref     = numel(hspace.ref);
end

C         = BezierExtractionOperator(basis,hspace);
ML        = MultiLevelOperator(hspace);
[B,shape] = BernsteinFunctionBasis(basis,hspace); 

hfuncs    = ML{hspace.level} * C{hspace.level} * B;
 
vals      = hpoints{nref}' * hfuncs;

if size(hpoints{nref},2) == 1
    X = reshape(vals(1,:),shape);
    Y = NaN;
    Z = NaN;
    
elseif size(hpoints{nref},2) == 2
    X = reshape(vals(1,:),shape);
    Y = reshape(vals(2,:),shape);
    Z = NaN;
    
elseif size(hpoints{nref},2) == 3
    X = reshape(vals(1,:),shape);
    Y = reshape(vals(2,:),shape);
    Z = reshape(vals(3,:),shape);
end
    
end