function [ML] = MultiLevelOperator(hspace)
% This function makes the multi-level operators for up to the highest
% refinement level for the hierarchical space defined in hspace. 
%
% Input:
%   hspace - data structure containing information on hierarchical basis
%
% Output:
%   ML - multi-level operators for each basis up to
%        hspace.level stored as cells for each level
%
%
% Created by Pieter van Zuijlen

assert(hspace.truncated == 0 | hspace.truncated == 1,'truncated input should be 0 for normal or 1 truncated');

R = hspace.R;

% get last refinement
nref = numel(hspace.ref);

% Build truncated refinement operator
if hspace.truncated == 1
    R_trunc{1}{1} = R{1}{1};
    for L = 1:hspace.level-1
        if any(hspace.ref{nref}.passive{L+1})
            J = diag(spones(hspace.ref{nref}.passive{L+1}));
            R_trunc{L}{L+1}   = R{L}{L+1} * J;
            R_trunc{L+1}{L+1} = R{L+1}{L+1};
        else
            % in this case there are no passive functions on this level,
            % and the regular R is taken.
            R_trunc{L}{L+1}   = R{L}{L+1};
            R_trunc{L+1}{L+1} = R{L+1}{L+1};
        end
    end   
    for L2 = 3:hspace.level
        for L1 = 1:L2-2
            R_trunc{L1}{L2} = R_trunc{L1}{L2-1} * R_trunc{L2-1}{L2};
        end
    end
    R = R_trunc;
end

for L = 1:hspace.level

    ML{L} = sparse([]);
    
    for K = 1:L    
        M = R{K}{L}(find(hspace.ref{nref}.active{K}),:);
        ML{L} = [ML{L}; M];
    end

end

end