function [hmesh] = MakeHmesh(basis,hspace)
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

elems    = basis.nelem/2;   %it will be multiplied with 2 in the for loop for lvl 1 as well
degree   = basis.degree;
nref     = numel(hspace.ref);

%% 2 dimensions

for L = 1:hspace.level
    
    elems = elems*2;
    funcs = elems+degree;
    
    hmesh{L} = zeros(1,prod(elems));
       
    elems_u = sort([ones(1,degree),1:elems(1),ones(1,degree)*elems(1)]);
    elems_v = sort([ones(1,degree),1:elems(2),ones(1,degree)*elems(2)]);
    
    elem_ind{L} = [];
    
    for i = 1:length(hspace.ref{nref}.deactivated{L})
        if hspace.ref{nref}.deactivated{L}(i) == 1
            
            iv = ceil(i/funcs(1));          %get funtion index in v direction
            iu = i - (iv-1)*funcs(1);       %get funtion index in u direction    
            
            ev = unique(elems_v(iv:iv+degree)); %get element indices in v direction
            eu = unique(elems_u(iu:iu+degree)); %get element indices in u direction
            
            for e = ev
                elem_ind{L} = [elem_ind{L} (eu+(e-1)*elems(1))];
            end
        end
    end
    hmesh{L}(elem_ind{L}) = 1;
    
end