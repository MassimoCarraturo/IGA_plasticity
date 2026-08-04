% BASISFUN_MULTI: Computes the multivariate tensor-product B-splines of
% degrees pl and with knot vectors u_knotl at the points contained in uu.
%
% INPUT
% pl:       vector of degrees in each direction
% uu:       matrix containing, in the j-th column, the j-th coordinate of the evaluation points
% u_knotl:  cell array containing the knot vectors
%
% OUTPUT
% bvalues:  values of the basis functions at the points (one column for each point)
%
% Author: Cesare Bracco
%
% Assembled from triplets in a single sparse() call.  The previous version
% built the result one column at a time (bvalues(:,h) = ...), which
% reallocated and copied the whole sparse array on every one of the npoints
% assignments; with the level-4 sphere (npoints = 155034, 11560 rows) that
% dominated the entire QI projection, costing minutes per call against the
% ~0.7 s the assembly actually needs.  Only the (p+1)^dim locally non-zero
% entries per point are formed, instead of a full-length kron per point per
% direction.  The result is bitwise identical to the previous implementation:
% the same two factors are multiplied per entry (IEEE multiplication is
% commutative) and sparse() drops structural zeros exactly as sparse(kron(...))
% did.  Verified with isequal() in 2-D and 3-D, on C^{p-1} and C^0 (interior
% knots of multiplicity p) spaces, with mixed degrees, non-uniform knots, and
% evaluation points on knots and at the parametric boundary.

function [bvalues] = basisfun_multi(pl, uu, u_knotl)

dim     = length(pl);
npoints = size(uu,1);

nfun  = zeros(1,dim);
spans = cell(1,dim);
Bd    = cell(1,dim);
for j = 1:dim
    nfun(j)  = length(u_knotl{j}) - pl(j) - 1;   % dimension in direction j
    sj       = findspan(nfun(j)-1, pl(j), uu(:,j).', u_knotl{j});
    spans{j} = sj(:);
    Bd{j}    = basisfun(sj, uu(:,j).', pl(j), u_knotl{j});   % npoints x (pl(j)+1)
end

% Tensor product over directions, keeping only the locally non-zero entries.
% Direction 1 varies fastest, matching kron(basisjj_j, basis_previous).
rows   = spans{1} - pl(1) + (0:pl(1));           % npoints x (pl(1)+1), 0-based
vals   = Bd{1};
stride = nfun(1);
for j = 2:dim
    base_j  = spans{j} - pl(j);
    nc      = size(rows,2);
    newRows = zeros(npoints, nc*(pl(j)+1));
    newVals = zeros(npoints, nc*(pl(j)+1));
    for b = 0:pl(j)
        idx            = b*nc + (1:nc);
        newRows(:,idx) = rows + (base_j + b) * stride;
        newVals(:,idx) = Bd{j}(:,b+1) .* vals;
    end
    rows   = newRows;
    vals   = newVals;
    stride = stride * nfun(j);
end

cols    = repmat((1:npoints).', 1, size(rows,2));
bvalues = sparse(rows(:)+1, cols(:), vals(:), prod(nfun), npoints);

end
