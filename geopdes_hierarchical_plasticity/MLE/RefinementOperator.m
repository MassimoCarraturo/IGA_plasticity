function [R] = RefinementOperator(basis, hspace)



assert(hspace.truncated == 0 | hspace.truncated == 1,'truncated input should be 0 for normal or 1 truncated');

% Define 1-d or n-d case

if iscell(basis.knots) == 0
    dim = 1;
    univar = 1;
else
    dim = numel(basis.knots);
    univar = 0;
end


% Loop over dimensions
for d = 1:dim

    levels  = hspace.nlevels;
    
    if univar
        knot{1} = basis.knots;
    else
        knot{1} = basis.knots{d};
    end
    
    p       = basis.order(d)-1;
    n       = length(knot{1})-p-1;
    % Make the first level to first level operator
    Rdim{d}{1}{1} = eye(n,n);


    for L = 1:levels

        %Define refined knots and new knotvector
        [~,~,newknots] = kntrefine (knot{L}, 1, p, 1);
        knot{L+1} = sort([knot{L} newknots]);

        n = length(knot{L}) - p - 1;
        m = length(newknots);

        % Construct the zero degree global insertion operator
        T{1} = zeros(n+m,n);
        for i = 1:n+m
            for j = 1:n
                if knot{L+1}(i) >= knot{L}(j) && knot{L+1}(i) < knot{L}(j+1)
                    T{1}(i,j) = 1;
                end
            end
        end

        % Construct the other degree global insertion operators
        for q = 1:p
            T{q+1} = zeros(n+m,n);

            % Add extra zero column to prevent index error
            T{q} = [T{q} zeros(n+m,1)];
            for i = 1:n+m
                for j = 1:n
                    % Split up left and right side for NaN check
                    left  = (knot{L+1}(i+q)-knot{L}(j))/(knot{L}(j+q)-knot{L}(j))*T{q}(i,j);
                    right = (knot{L}(j+q+1)-knot{L+1}(i+q))/(knot{L}(j+q+1)-knot{L}(j+1))*T{q}(i,j+1);
                    if isnan(left) == 1
                        left = 0;
                    elseif isnan(right) == 1
                        right = 0;
                    end
                    T{q+1}(i,j) = left + right;
                end
            end
        end

        % The last global insertion operator is the Refinement Operator of that
        % level
        Rdim{d}{L}{L+1}   = T{q+1}';
        % Define L+1th level to L+1th level operator
        Rdim{d}{L+1}{L+1} = eye(n+m,n+m); 
    end

    for L2 = 3:levels+1
        for L1 = 1:L2-2
            Rdim{d}{L1}{L2} = Rdim{d}{L1}{L2-1} * Rdim{d}{L1+1}{L2};
        end
    end
end

R = Rdim{1};
if univar == 0
    for d = 2:dim
        for L2 = 1:levels+1
            for L1 = 1:L2
                R{L1}{L2} = kron(Rdim{d}{L1}{L2},R{L1}{L2});
            end
        end
    end
end


% Build truncated refinement operator
if hspace.truncated == 1
    R_trunc{1}{1} = R{1}{1};
    for L = 1:levels
        R_trunc{L}{L+1}   = R{L}{L+1} * diag(hspace.passive{L+1});
        R_trunc{L+1}{L+1} = R{L+1}{L+1};
    end   
    for L2 = 3:levels+1
        for L1 = 1:L2-2
            R_trunc{L1}{L2} = R_trunc{L1}{L2-1} * R_trunc{L1+1}{L2};
        end
    end
    R = R_trunc;
end

end