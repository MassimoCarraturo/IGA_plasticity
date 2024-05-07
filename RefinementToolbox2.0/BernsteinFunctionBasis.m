function [B,shape] = BernsteinFunctionBasis(basis,hspace,level)
% Creation of the Bernstein polynomial basis of input level
%
% Input:
%   basis  - data structure containing information on the parametric basis
%   geom   - data structure containing information on the geometric basis
%   hspace - data structure containing information on the hierarchcal basis
%
%   (optional) level  - level of Bernstein basis. Without this value, the
%   level of hspace is taken.
%
% Output:
%   B      - set of functions of Bernstein basis.
%   shape  - amount of values for all space dimensions. Is used for
%   reshape.
%
%
% Created by Pieter van Zuijlen

if ~exist('level','var')
    level = hspace.level;
end

degree = basis.degree;
nelem  = basis.nelem*2^(level-1);
n      = round(basis.resol./nelem);
dim    = basis.dim;

assert(dim < 3, 'a dimension higher than 2 is not possible');
assert(length(n) == length(nelem), 'dimension error')

B = 1;

for d = 1:dim

    Belem = zeros(n(d),degree+1);
    i = 0;
    for xi = linspace(-1,1,n(d)+1)
        i      = i + 1;
        Belem(i,:) = brnstpoly(xi,degree);
    end

     Bend         = Belem(end,:);
     Belem(end,:) = [];

    nbpoly = nelem(d) * (degree+1);
    Bd  = sparse(nelem(d)*n(d)+1,nbpoly);

    for e = 1:nelem(d)
        ielem = (e-1)*(degree+1)+1:e*(degree+1);
        ixi   = (e-1)*n(d)+1:e*n(d);
        Bd(ixi,ielem) = Belem;
    end

    Bd(end,end-degree:end) = Bend;
    
    Bd = Bd';

    B = kron(Bd,B);

end

shape = n.*nelem+1;

if dim == 1
    shape = [1,shape];
end





% if ~exist('level','var')
%     level = hspace.level;
% end
% 
% degree = basis.degree;
% nelem  = basis.nelem*2^(level-1);
% n      = round(geom.res./nelem);
% dim    = basis.dim;
% 
% assert(dim < 3, 'a dimension higher than 2 is not possible');
% assert(length(n) == length(nelem), 'dimension error')
% 
% if dim == 1;
%     
%     Belem = zeros(n,degree+1);
%     i = 0;
%     for xi = linspace(-1,1,n+1)
%         i      = i + 1;
%         Belem(i,:) = brnstpoly(xi,degree);
%     end
%     
%      Bend         = Belem(end,:);
%      Belem(end,:) = [];
%    
%     nbpoly = nelem * (degree+1);
%     B  = sparse(nelem*n+1,nbpoly);
% 
%     for e = 1:nelem
%         ielem = (e-1)*(degree+1)+1:e*(degree+1);
%         ixi   = (e-1)*n+1:e*n;
%         B(ixi,ielem) = Belem;
%     end
%     
%      B(end,end-degree:end) = Bend;
%     
%     B = B';
%     
%     shape = size(B,2);
% 
% elseif dim == 2;
%    
%     for d = 1:dim
%         Bd{d} = zeros(n(d),(degree+1));
%         i = 0;
%         for xi = linspace(-1,1,n(d))
%             i      = i + 1;
%             Bd{d}(i,:) = brnstpoly(xi,degree);
%         end 
%     end
% 
%     Belem = zeros(prod(n),(degree+1)^dim);
%     ind_mat = reshape(1:(degree+1)^dim,[degree+1,degree+1]);
% 
%     for i = 1:degree+1
%         for j = 1:degree+1
% 
%             ind = ind_mat(i,j);
%             Belem(:,ind) = kron(Bd{1}(:,i),Bd{2}(:,j));
% 
%         end
%     end
%     
%     
%     nfunc  = size(Belem,2);
%     order  = nfunc ^ (1/dim);
%     degree = order-1;
% 
%     xi = zeros(nelem.*[n(1),n(2)]);
% 
%     B = sparse(n(1)*n(2)*prod(nelem),prod(nelem)*order^2);
%     count = 0
%     
%     for e2 = 1:nelem(2)
%         for i = 1:order
%             for e1 = 1:nelem(1)
%                 for j = 1:order
%                     count = count + 1;
%                     ifun  = j + order * (i-1);
%                     iu    = (e1-1)*n(1)+1:e1*n(1);
%                     iv    = (e2-1)*n(2)+1:e2*n(2);
% 
%                     func        = xi;
%                     func(iv,iu) = reshape(Belem(:,ifun),[n(2),n(1)]);
%                     B(:,count)  = reshape(func,[n(1)*n(2)*prod(nelem),1]);
%                 end
%             end
%         end
%     end
%     
%     shape = nelem.*n;
%     
% elseif dim == 3
%     
%     for d = 1:dim
%         Bd{d} = zeros(n(d),(degree+1));
%         i = 0;
%         for xi = linspace(-1,1,n(d))
%             i      = i + 1;
%             Bd{d}(i,:) = brnstpoly(xi,degree);
%         end 
%     end
%     
%     B = zeros(prod(n),(degree+1)^dim);
%     ind_mat = reshape(1:(degree+1)^dim,[degree+1,degree+1,degree+1]);
%     
%     for i = 1:degree+1
%         for j = 1:degree+1
%             for k = 1:degree+1
%                 ind = ind_mat(i,j,k);
% 
%                 B(:,ind) = kron(kron(Bd{3}(:,i),Bd{2}(:,j)),Bd{1}(:,k));
%             end
%         end
%     end
%     
%     nfunc  = size(Belem,2);
%     order  = nfunc ^ (1/dim);
%     degree = order-1;
% 
%     xi = zeros(nelem.*[n(1),n(2),n(3)]);
% 
%     B = sparse(n(1)*n(2)*n(3)*prod(nelem),prod(nelem)*order^2);
%     count = 0;
% 
%     for e3 = 1:nelem(3)
%         for k = 1:order
%             for e2 = 1:nelem(2)
%                 for i = 1:order
%                     for e1 = 1:nelem(1)
%                         for j = 1:order
%                             count = count + 1;
%                             ifun = j + order*(i-1)+order^2*(k-1);
%                             iu   = (e1-1)*n(1)+1:e1*n(1);
%                             iv   = (e2-1)*n(2)+1:e2*n(2);
%                             iw   = (e3-1)*n(3)+1:e3*n(3);
% 
%                             func = xi;
%                             func(iu,iv,iw) = reshape(Belem(:,ifun),[n(1),n(2),n(3)]);
%                             B(:,count) = reshape(func,[n(1)*n(2)*n(3)*prod(nelem),1]);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end    



end