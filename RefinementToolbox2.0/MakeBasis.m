function [basis] = MakeBasis(nurbs,resolution)
% Creation of the basis data struct with the NURBS toolbox data structure. 
% Basis contains information on the parametrical space. If you are not 
% using the NURBS Toolbox, use function CreateGeom.m instead.
%
% Input:
%   NURBS - Data struct created in the NURBS Toolbox
%
% Output:
%   basis  - Parametrical data structure
%
%   basis.dim    - dimension of the parametrical space
%   basis.knot   - knot vectors for each direction stored as cells
%   basis.degree - polynomial degree
%   basis.nfunc  - amount of basis functions spanning the parametric space
%   basis.nelem  - list of amount of elements in each direction
%   basis.resol  - amount of points used to discretize the geometry
%
% Created by Pieter van Zuijlen

if iscell(nurbs.knots) == 0;
    basis.dim        = 1;
    basis.knot{1}    = nurbs.knots;
    basis.degree     = nurbs.order - 1;
    basis.nfunc      = length(basis.knot{1}) - basis.degree - 1;
    basis.nelem      = basis.nfunc - basis.degree;
    basis.resol      = resolution(1);
else
    basis.dim     = numel(nurbs.knots);
    basis.degree  = nurbs.order(1)-1;
    if length(resolution) == basis.dim
    	basis.resol = resolution;
    else
        basis.resol = ones(1,basis.dim)*resolution(1);
    end
    for d = 1:basis.dim
        basis.knot{d} = nurbs.knots{d};
        if d >= 2
            assert(nurbs.order(d) == nurbs.order(d-1),'The function order should be the same in each direction');
        end
        basis.nfunc(d)   = length(basis.knot{d}) - basis.degree - 1;
        basis.nelem(d)   = basis.nfunc(d) - basis.degree;
    end
end


end