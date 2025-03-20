function [MLC] = MultiLevelExtractionOperator(basis, hspace, MultLvl)

% Define if 1-d or n-d case

if iscell(basis.knots) == 0
    dim = 1;
    univar = 1;
else
    dim = numel(basis.knots);
    univar = 0;
end

levels  = hspace.level;
order   = basis.order(1);
p       = order-1;


for L = 1:levels+1
    
    C = 1;
    
    for d = 1:dim
        
        if univar
            knot{d}{1} = basis.knots;
        else
            knot{d}{1} = basis.knots{d};
        end
    
        n       = length(knot{d}{1})-p-1;

        [~,~,newknots] = kntrefine (knot{d}{L}, 1, p, 1);
        knot{d}{L+1} = sort([knot{d}{L} newknots]);

        n      = length(knot{d}{L})-p-1;
        e      = n - p;
        nb     = e*order;

        ML = MultLvl{L};

        Ce = bzrextr(knot{d}{L},p);

        Cdim{d} = sparse(n,nb);

        for i = 1:n-p
            Cdim{d}(i:i+p,(i-1)*order+1:i*order) = Ce(:,:,i);
        end
        
        C = kron(C,Cdim{d});
        
    end
    
    MLC{L} = ML * C;

end

end

% levels  = hspace.level;
% knot{1} = basis.knots;
% order   = basis.order(1);
% p       = order-1;
% n       = length(knot{1})-p-1;
% 
% for L = 1:levels+1
% 
%     [~,~,newknots] = kntrefine (knot{L}, 1, p, 1);
%     knot{L+1} = sort([knot{L} newknots]);
%     
%     n      = length(knot{L})-p-1;
%     e      = n - p;
%     nb     = e*order;
% 
%     ML = MultLvl{L};
% 
%     Ce = bzrextr(knot{L},p);
% 
%     C = sparse(n,nb);
% 
%     for i = 1:n-p
%         C(i:i+p,(i-1)*order+1:i*order) = Ce(:,:,i);
%     end
%     
%     MLC{L} = ML * C;

