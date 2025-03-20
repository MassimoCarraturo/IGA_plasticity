function [ML] = MultiLevelOperator(hspace,R)

levels = hspace.level;

for L = 1:levels+1

    ML{L} = [];
    
    for K = 1:L
        
        % Make row inclusion matrix
        J = diag(hspace.active{K});
        J( ~any(J,2), : ) = []; % set all zero rows to empty 
               
        M = J * R{K}{L};
        ML{L} = [ML{L}; M];
    end

end

end