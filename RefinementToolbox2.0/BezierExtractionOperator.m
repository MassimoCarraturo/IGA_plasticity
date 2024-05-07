function [C] = BezierExtractionOperator(basis, hspace)
% This function makes Bezier Extraction operators for each basis up to
% the refinement level. The function uses kntrefine.m and bzrextr of the
% NURBS toolbox
%
% Input:
%   basis  - data structure containing information on the parametric bais
%   hspace - data structure containing information on hierarchical basis
%
% Output:
%   C  - struct with Bezier extraction operators for each basis up to
%         hspace.level
%
%
% Created by Pieter van Zuijlen

p       = basis.degree;
order   = p+1;


for L = 1:hspace.level
    
    C{L} = 1;
    
    for d = 1:basis.dim
        
        knot{d}{1} = basis.knot{d};

        n       = length(knot{d}{1})-p-1;

        [~,~,newknots] = kntrefine (knot{d}{L}, 1, p, p-1);
        knot{d}{L+1} = sort([knot{d}{L} newknots]);
        
        
        n      = length(knot{d}{L})-p-1;
        e      = n - p;
        nb     = e*order;

        Ce = bzrextr(knot{d}{L},p);

        Cdim{d} = sparse(n,nb);

        for i = 1:n-p
            Cdim{d}(i:i+p,(i-1)*order+1:i*order) = Ce(:,:,i);
        end
        
        C{L} = kron(Cdim{d},C{L});
        
    end

end

end


