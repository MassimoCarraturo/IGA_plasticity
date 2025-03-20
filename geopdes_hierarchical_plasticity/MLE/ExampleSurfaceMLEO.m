clear all, close all

load('ExamSurf')

srf = ExamSurf.srf;

degree = srf.order(1)-1;
truncated = 1;

nu = 10;
nv = 10;


%% Defining multivariate hspace

hspace.level          = 2;

hspace.truncated      = 0;

hspace.active{1}      = [0 0 0 0 0 0 0 0 0];
hspace.active{2}      = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
hspace.active{3}      = [1 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

hspace.deactivated{1} = [1 1 1 1 1 1 1 1 1];
hspace.deactivated{2} = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
hspace.deactivated{3} = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

hspace.passive{1}     = [1 1 1 1 1 1 1 1 1];
hspace.passive{2}     = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
hspace.passive{2}     = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
hspace.passive{3}     = [0 0 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

%% origional surface

[X,Y,Z] = MakeSurf(nu*4,nv*4,srf);

cpts = reshape(srf.coefs,[4,prod(srf.number)]);

figure
subplot(2,2,1)
h1 = surf(X,Y,Z)
title('Original surface')
hold on
scatter3(cpts(1,:)',cpts(2,:)',cpts(3,:)','r','filled')
axis([-1,1,-1,1,0,1])
view(-20,20)
set(h1, 'edgecolor','none')

subplot(2,2,2)
h2 = surf(X,Y,Z)
title('Hierarchical structure mesh')
hold on
plot3(X(:,1),Y(:,1),Z(:,1),'LineWidth',2,'Color','k')
plot3(X(:,end),Y(:,end),Z(:,end),'LineWidth',2,'Color','k')
plot3(X(1,:),Y(1,:),Z(1,:),'LineWidth',2,'Color','k')
plot3(X(end,:),Y(end,:),Z(end,:),'LineWidth',2,'Color','k')
plot3(X(:,nu*2),Y(:,nu*2),Z(:,nu*2),'LineWidth',2,'Color','k')
plot3(X(nv*2,:),Y(nv*2,:),Z(nv*2,:),'LineWidth',2,'Color','k')
plot3(X(1:nv*2,nu),Y(1:nv*2,nu),Z(1:nv*2,nu),'LineWidth',2,'Color','k')
plot3(X(nv,1:nu*2),Y(nv,1:nu*2),Z(nv,1:nu*2),'LineWidth',2,'Color','k')
axis([-1,1,-1,1,0,1])
view(-20,20)
set(h2, 'edgecolor','none')


%% 

R = RefinementOperator(srf,hspace);

ML  = MultiLevelOperator(hspace,R);

MLC = MultiLevelExtractionOperator(srf,hspace,ML);

K = RefinementMatrix(hspace,R);

%% construct bernstein polynomial basis

Belem = brnstbasis([nu,nv],degree);
elems = [4 4];

B = constbspc(Belem, elems, nu, nv);

cpts = reshape(srf.coefs,[4,prod(srf.number)]);
hcpts= K * cpts';

bcpts= MLC{3}' * hcpts; 

p = bcpts' * B';

X = reshape(p(1,:),elems.*[nu,nv]);
Y = reshape(p(2,:),elems.*[nu,nv]);
Z = reshape(p(3,:),elems.*[nu,nv]);

subplot(2,2,3)
h3 = surf(X,Y,Z)
title('Refined surface')
hold on
scatter3(hcpts(:,1),hcpts(:,2),hcpts(:,3),'r','filled')
axis([-1,1,-1,1,0,1])
view(-20,20)
set(h3, 'edgecolor','none')

HB = MLC{3} * B';



p = hcpts' * HB;

X = reshape(p(1,:),elems.*[nu,nv]);
Y = reshape(p(2,:),elems.*[nu,nv]);
Z = reshape(p(3,:),elems.*[nu,nv]);

subplot(2,2,4)
h4 = surf(X,Y,Z)
title('Refined surface with Bernstein basis')
hold on
scatter3(bcpts(:,1),bcpts(:,2),bcpts(:,3),'r','filled')
axis([-1,1,-1,1,0,1])
view(-20,20)
set(h4, 'edgecolor','none')

% 
% elem = srf.number-srf.order+1;
% 
% levels = hspace.level;
% 
% elem_ref = elem * 2^levels;
% 
% bcpts = reshape(bcpts',[4,12,12])
% 
% bknot{1} = sort([linspace(0,4,5) linspace(0,4,5) linspace(0,4,5)]);
% bknot{2} = sort([linspace(0,4,5) linspace(0,4,5) linspace(0,4,5)]);
% 
% bsrf = nrbmak(bcpts,bknot);
% 
% NV = 100;
% NU = 100;
% 
% [SURF,IND] = nrbbasisfun({linspace(0.0,1.0,NV), linspace(0.0,1.0,NU)}, bsrf);
% 
% INDligne=[1:NU*NV]'*ones(1,bsrf.order(1)*bsrf.order(2));
% SURF1=full(sparse(INDligne,IND,SURF));
% 
% elems = [4,4];
% 
% Belem = brnstbasis([nu,nv],degree);
% 
% B1 = [];
% B2 = [];
% B3 = [];
% B = [];
% 
% for i = 1:4
%     B1 = [B1 SURF(:,1:3)];
%     B2 = [B2 SURF(:,4:6)];
%     B3 = [B3 SURF(:,7:9)];
% end
% 
% for i = 1:4
%     B = [B B1 B2 B3];
% end
% 
% HB = MLC{3} * B';
% 
% coefs = reshape(srf.coefs,[4,9]);
% 
% hcpts = K * coefs';
% 
% 
% Pij=reshape(bsrf.coefs,4,[])';
% p = reshape(B*Pij,[NV,NU,4]);
% 
% X(:,:) = p(:,:,1);
% Y(:,:) = p(:,:,2);
% Z(:,:) = p(:,:,3);

