function [K] = RefinementOperator(hspace,refnbr)
% Creation of the Hierarchial level operator, for truncated or normal
% case. The operator is defined from the entered refinement number to the
% latest refinement. If no entry is given, the unrefined basis is taken.
%
% Input:
%   hspace - data structure containing information on hierarchical basis
%   (optional) refnbr - map the values from the hierarchical space of this
%   refinement to the latest refinement.
%
% Output:
%   K  - Hierarchical level operator. The matrix expressing control points
%        of the hierarchical level as a linear combination of control 
%        points of level 1.
%
%
% Created by Pieter van Zuijlen

if ~exist('refnbr','var')
    refnbr = 1;
end

nref = numel(hspace.ref);

old = hspace.ref{refnbr};
new = hspace.ref{nref};

if hspace.truncated == 0
    
    % Get the matrix for the first level. This is an identity matrix for a
    % unrefined mesh   
    J{1} = speye(length(old.active{1}));
    J{1} = J{1}(find(new.active{1}+new.deactivated{1}),find(old.active{1}));

    K{1} = J{1};
          
    for L1 = 1:hspace.level-1
               
        iAct = [];
        for i = 1:L1-1
            iAct = [iAct hspace.ref{nref}.active{i}(find(hspace.ref{nref}.active{i}))];
        end
        iAct = [iAct hspace.ref{nref}.active{L1}(find(hspace.ref{nref}.deactivated{L1}+hspace.ref{nref}.active{L1}))];


        iD     = find(new.deactivated{L1});
        iAD    = find(new.active{L1+1} + new.deactivated{L1+1});

        K_A = K{L1}(find(iAct),:);
        K_D = K{L1}(find(~iAct),:);

        Rsub   = hspace.R{L1}{L1+1}(iD,iAD)';
        
        % Check if level L already existed in the old hspace
        if numel(old.active) < L1+1
            K{L1+1} = [K_A; Rsub*K_D];
        else        
            J{L1+1} = eye(length(old.active{L1+1}));
            J{L1+1} = J{L1+1}(find(new.active{L1+1}+new.deactivated{L1+1}),find(old.active{L1+1}));
            O      = zeros(size(K_A,1),size(J{L1+1},2));
            K{L1+1} = [K_A, O; Rsub*K_D, J{L1+1}];
        end

    end

    K     = K{hspace.level};

else
      
    J{1} = speye(length(old.active{1}));
    J{1} = J{1}(:,find(old.active{1})); 

    Ksub = J{1};

    K    = Ksub(find(new.active{1}),:);

    for L = 1:hspace.level-1

        Rsub = hspace.R{L}{L+1}';

        if numel(old.active) >= L+1
            J{L+1} = speye(length(old.active{L+1}));
            J{L+1} = J{L+1}(:,find(old.active{L+1}));
            
%           Erase rows in R for information that is already in the previous
%           state
            Rsub(find(old.active{L+1}),:) = zeros(length(find(old.active{L+1})),length(old.active{L}));

        else 
            J{L+1} = [];
        end

        Ksub = [Rsub * Ksub , J{L+1}];
        O    = zeros(size(K,1),size(Ksub,2)-size(K,2));
        K    = [K , O ; Ksub(find(new.active{L+1}),:)];

    end
    
end





end