function [hspace] = RefineHspace(basis, hspace, refFunc)
% Update the hspace after refining the selected functions. 
%
% Input:
%   basis     - data structure containing information on parametric basis
%   hspace    - data structure containing information on hierarchical basis
%   refFunc   - list of function indices of functions to be refined
%
% Output:
%   hspace    - Updated with the refined function structure
%
%
% Created by Pieter van Zuijlen

for fun = refFunc
    assert(fun <= hspace.nfunc,'refFunc is higher than the amount of functions in hspace')
end

nFunc = length(refFunc);
iSpac = [];

% Store the old hspace data

nref = numel(hspace.ref);

hspace.ref{nref+1}.active      = hspace.ref{nref}.active;
hspace.ref{nref+1}.passive     = hspace.ref{nref}.passive;
hspace.ref{nref+1}.deactivated = hspace.ref{nref}.deactivated;

for iFunc = refFunc
    n = 0;
    for L = 1:hspace.level
        iAct = find(hspace.ref{nref}.active{L});
        n = n + length(iAct);
        if n >= iFunc
            iLvl  = iAct(iFunc-n+length(iAct)); 
            iSpac = [iSpac; iLvl, L];
            break
        end
    end
end

for i = 1:nFunc
    
    iLvl = iSpac(i,1);
    L    = iSpac(i,2);
    
    hspace.ref{nref+1}.active{L}(iLvl) = 0;
    hspace.ref{nref+1}.deactivated{L}(iLvl) = 1;

    %update R and hspace if a new level of refinement is reached
    if hspace.level == L
        hspace.level = L+1;
        tic
        hspace.R = LevelOperator(basis,hspace);
        toc
        
        for d = 1:basis.dim
            e(d) = length(basis.knot{d})-2*basis.degree-1;
        end
        nNewFunc = prod(e*2^L + basis.degree);%amount of new functions
        hspace.ref{nref+1}.active{L+1}      = zeros(1,nNewFunc);
        hspace.ref{nref+1}.deactivated{L+1} = zeros(1,nNewFunc);
        hspace.ref{nref+1}.passive{L+1}     = ones(1,nNewFunc);
        
    end

    iLvlPlus = find(hspace.R{L}{L+1}(iLvl,:));
    for index = iLvlPlus
        % If the function has already been deactivated, it can not become
        % active again.
        if hspace.ref{nref+1}.deactivated{L+1}(index) == 0
            hspace.ref{nref+1}.active{L+1}(index)  = 1;
        end
    end
    hspace.ref{nref+1}.passive{L+1}(iLvlPlus) = 0;   
    
    % Check if not all passive function passing information to the active
    % functions are gone. (for truncated case)
    if hspace.truncated == 1
        for index = find(hspace.ref{nref+1}.active{L})
            iLvlPlus = find(hspace.R{L}{L+1}(index,:));
            if sum(hspace.ref{nref+1}.passive{L+1}(iLvlPlus)) == 0
                hspace.ref{nref+1}.active{L}(index) = 0;
                hspace.ref{nref+1}.deactivated{L}(index) = 1;
            end
        end
    end
    
end

hspace.nfunc = 0;

for L = 1:hspace.level
    hspace.nfunc = hspace.nfunc + length(find(hspace.ref{nref+1}.active{L}));
end

end