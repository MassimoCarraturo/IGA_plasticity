function [K] = RefinementMatrix(hspace, R)

levels  = hspace.level;

% Get the matrix for the first level. This is an identity matrix for a
% unrefined mesh
J{1} = diag(hspace.active{1} + hspace.deactivated{1});
J{1}( ~any(J{1},2) , : ) = [];
J{1}( : , ~any(J{1},2) ) = [];

K{1} = J{1};

if hspace.truncated == 0
    for L = 1:levels

        act_hsp = [];

        for i = 1:L-1
            act_hsp = [act_hsp hspace.active{i}(find(hspace.active{i}))];
        end
        act_hsp = [act_hsp hspace.active{L}(find(hspace.deactivated{L}+hspace.active{L}))];

        deact_lvl = find(hspace.deactivated{L});
        act_and_deact_lvl = find(hspace.active{L+1} + hspace.deactivated{L+1});



        % Is non empty for applications not beginning at unrefined mesh
        J{L+1} = [];

        K_A = K{L}(find(act_hsp),:);
        K_D = K{L}(find(~act_hsp),:);

        Rsub = R{L}{L+1}(deact_lvl,act_and_deact_lvl)';

        K{L+1} = [K_A; Rsub*K_D];

    end
    
    K = K{levels+1};
    
else
    for L = 1:levels

        act_hsp = [];
        for i = 1:L-1
            act_hsp = [act_hsp hspace.active{i}(find(hspace.active{i}))];
        end
        act_hsp = [act_hsp hspace.active{L}]

        % Is non empty for applications not beginning at unrefined mesh
        J{L+1} = [];
        
        K_A = K{L}(find(act_hsp),:);
        
        if L == 1
            K_{L} = K{1};
        else
            K_{L} = R{L-1}{L}' * K_{L-1};
        end
        
        size(K_A)
        size(R{L}{L+1}')
        size(K_{L})
        
        K{L+1} = [K_A; R{L}{L+1}'*K_{L}];

    end
    
%     K = K{levels+1}(find(act_hsp),:);
    
end


end