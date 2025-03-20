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

function [bvalues] = basisfunder_multi(pl, uu, u_knotl)

dim=length(pl);
npoints=size(uu,1);
basis=cell(1,npoints);
for h=1:npoints
    basis{h}=1;
end

for j=1:dim
    dofj=length(u_knotl{j})-pl(j)-1;
    for h=1:npoints
        ii=findspan(dofj-1,pl(j),uu(h,j),u_knotl{j});
        basisjj(:,h)=zeros(dofj,1);
        basisjj(ii-pl(j)+1:ii+1,h)=basisfun (ii, uu(h,j), pl(j), u_knotl{j});
        basis{h}=sparse(kron(basisjj(:,h),basis{h}));
    end
end

for h=1:npoints
    aux=basis{h};
    bvalues(:,h)=aux(:);
end

end