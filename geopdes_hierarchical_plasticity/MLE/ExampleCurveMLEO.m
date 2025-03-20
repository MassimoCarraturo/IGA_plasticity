close all, clear all

%% Defining space basis

closed    = [0 1 2 3 4];
degree    = 2;
nu        = 1000;

knot   = kntbrkdegreg (closed, degree);
n      = length(knot) - degree - 1;
Py     = [1 4 2 4 6 3 2 4 5 6 7 3 2 4 5 ];
Py     = Py(1:n);
Px     = linspace(0,4,n);
e      = n - degree;

xi     = linspace(knot(1),knot(end),nu);

crv = nrbmak(Py,knot);

%% Defining hierarchical refined structure

hspace.level          = 2;

hspace.truncated      = 1;

hspace.active{1}      = [1 1 1 1 0 0];
hspace.active{2}      = [0 0 0 0 0 0 1 0 0 0];
hspace.active{3}      = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1];

hspace.deactivated{1} = [0 0 0 0 1 1];
hspace.deactivated{2} = [0 0 0 0 0 0 0 1 1 1];
hspace.deactivated{3} = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

hspace.passive{1}     = [0 0 0 0 0 0];
hspace.passive{2}     = [1 1 1 1 1 1 0 0 0 0];
hspace.passive{3}     = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0];

%% Multi level extraction operator 

R   = RefinementOperator(crv,hspace);

ML  = MultiLevelOperator(hspace,R);

MLC = MultiLevelExtractionOperator(crv,hspace,ML);


%% Control point calculation

hspace.truncated = 0;

%R   = RefinementOperator(crv,hspace);

K = RefinementMatrix(hspace,R);

hcptsy = K * Py';
hcptsx = K * Px';

Bcptsy = MLC{3}' * hcptsy;
Bcptsx = MLC{3}' * hcptsx;

%% Make Bernstein polynomial basis

[N I] = nrbbasisfun (xi, crv);

vals1y = Py(I) .* N;
valslx = Px(I) .* N;

bbzr1y = zeros(nu,1);
bbzr1x = zeros(nu,1);

for i = 1:crv.order(1)
    bbzr1y = bbzr1y+ vals1y(:,i);
    bbzr1x = bbzr1x+ valslx(:,i);
end





Bknot = sort([linspace(0,4,17),linspace(0,4,17),linspace(0,4,17)]);
Bcptsy = [Bcptsy zeros(48,2) ones(48,1)];

Bern = nrbmak(Bcptsy',Bknot);
[B I] = nrbbasisfun (xi, Bern);

vals2y = Bcptsy(I) .* B;
bbzr2y = zeros(nu,1);

vals2x = Bcptsx(I) .* B;
bbzr2x = zeros(nu,1);

for i = 1:crv.order(1)
    bbzr2y = bbzr2y+ vals2y(:,i);
    bbzr2x = bbzr2x+ vals2x(:,i);
end



%% Extraction from bernstein polynomials.

Bez = zeros(48,nu);

figure
for i = 1:48
    [ind1 ind2] = find((I) == i);
    ind2 = find((I) == i);
    Bez(i,ind1) = B(ind2)';
    plot(xi,Bez(i,:))
    hold on
end

HB = MLC{3} * Bez;

valsy = hcptsy' * HB;
valsx = hcptsx' * HB;

figure
subplot(2,3,1)
plot(bbzr1x,bbzr1y)
title('Unrefined curve')
hold on
scatter(Px,Py)
axis([0,4,0,6])

subplot(2,3,2)
plot(valsx,valsy)
title('Refined curve')
hold on
scatter(hcptsx,hcptsy)
axis([0,4,0,6])

subplot(2,3,3)
plot(bbzr2x,bbzr2y)
title('Refined curve with Bernstein basis')
hold on
scatter(Bcptsx,Bcptsy(:,1))
axis([0,4,0,6])

subplot(2,3,4)
plot(N,'b')
title('Original spline basis')

subplot(2,3,5)
plot(HB','b')
title('Refined spline basis')

subplot(2,3,6)
plot(B,'b')
title('Bernstein Polynomials')

% 
% figure
% plot(bbzr1x,bbzr1y)
% title('Refined curve')
% hold on
% scatter(Px,Py)
% 
% 
% figure
% plot(valsx,valsy)
% hold on
% scatter(hcptsx,hcptsy)
% 
% figure
% plot(bbzr2x,bbzr2y)
% title('Bezier constructed')
% hold on
% scatter(Bcptsx,Bcptsy(:,1))
% 
% figure
% plot(bbzr1x,bbzr1y)
% title('Unrefined')
% hold on
% scatter(Px,Py)
% 



%% Truncated

% figure
% subplot(1,2,1)
% plot(HB','b')
% title('Normal')
% 
% R_{1}{2} = R{1}{2}*diag(hspace.passive{2});
% R_{2}{3} = R{2}{3}*diag(hspace.passive{3});
% 
% R_{1}{3} = R_{1}{2}*R_{2}{3};
% 
% R_{1}{1} = R{1}{1};
% R_{2}{2} = R{2}{2};
% R_{3}{3} = R{3}{3};
% 
% ML  = MultiLevelOperator(hspace,R_);
% 
% MLC = MultiLevelExtractionOperator(crv,hspace,ML);
% 
% THB = MLC{3} * Bez; 
% 
% subplot(1,2,2)
% plot(THB','b')
% title('Truncated')

K3 =[1.0000         0         0         0         0         0         0 0 0 0 0;
         0    1.0000         0         0         0         0         0 0 0 0 0;
         0         0    1.0000         0         0         0         0 0 0 0 0;
         0         0         0    1.0000         0         0         0 0 0 0 0;
         0         0         0         0    1.0000         0         0 0 0 0 0;
    1.0000         0         0         0         0         0         0 0 0 0 0;
    0.7500    0.2500         0         0         0         0         0 0 0 0 0;
    0.3750    0.5625    0.0625         0         0         0         0 0 0 0 0;
    0.1250    0.6875    0.1875         0         0         0         0 0 0 0 0;
         0    0.6250    0.3750         0         0         0         0 0 0 0 0;
         0    0.3750    0.6250         0         0         0         0 0 0 0 0;
         0    0.1875    0.7500    0.0625         0         0         0 0 0 0 0;
         0    0.0625    0.7500    0.1875         0         0         0 0 0 0 0;
         0         0    0.6250    0.3750         0         0         0 0 0 0 0;
         0         0    0.3750    0.6250         0         0         0 0 0 0 0;
         0         0    0.1875    0.5625    0.2500         0         0 0 0 0 0;
         0         0    0.0625    0.1875    0.7500         0         0 0 0 0 0;
         0         0         0         0         0    1.0000         0 0 0 0 0;
         0         0         0         0         0         0    1.0000 0 0 0 0;
         0         0         0         0         0         0         0    1.0000         0         0         0;
         0         0         0         0         0         0         0         0    1.0000         0         0;
         0         0         0         0         0         0         0         0         0    1.0000         0;
         0         0         0         0         0         0         0      0         0         0    1.0000];


         