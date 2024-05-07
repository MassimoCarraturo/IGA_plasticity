function [R] = LevelOperator(basis, hspace)
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


% Loop over dimensions
for d = 1:basis.dim
    knot{1} = basis.knot{d};


    p       = basis.degree;
    n       = length(knot{1})-p-1;
    % Make the first level to first level operator
    Rdim{d}{1}{1} = speye(n,n);


    for L = 1:hspace.level-1

        %Define refined knots and new knotvector
        [~,~,newknots] = kntrefine (knot{L}, 1, p, p-1);
        knot{L+1} = sort([knot{L} newknots]);

        n = length(knot{L}) - p - 1;
        m = length(newknots);

        % Construct the zero degree global insertion operator
        T{1} = zeros(n+m,n);
        for i = 1:n+m
            for j = 1:n
                if knot{L+1}(i) >= knot{L}(j) & knot{L+1}(i) < knot{L}(j+1)
                    T{1}(i,j) = 1;
                end
            end
        end

        % Construct the other degree global insertion operators
        for q = 1:p
            T{q+1} = sparse(n+m,n);

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

    for L2 = 3:hspace.level
        for L1 = 1:L2-2
            Rdim{d}{L1}{L2} = Rdim{d}{L1}{L2-1} * Rdim{d}{L2-1}{L2};
        end
    end
end

R = Rdim{1};
if basis.dim ~= 1
    for d = 2:basis.dim
        for L2 = 1:hspace.level
            for L1 = 1:L2
                R{L1}{L2} = kron(Rdim{d}{L1}{L2},R{L1}{L2});
            end
        end
    end
end
hspace.R = R;    

end