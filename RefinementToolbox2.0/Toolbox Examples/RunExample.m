clear all, close all

%% Load data

load('Curve')

%% Settings

truncated = 1;
nu        = 256;
nv        = 256;
dim       = 2;

%% Initiation

basis   = MakeBasis(nurbs,[nu,nv]);
hspace  = MakeHspace(nurbs, truncated);
hpoints = MakeHpoints(nurbs,dim);

%% Construct unrefined geometry

[X0,Y0,~,hfuncs0,shp] = ConstructGeometry(basis,hspace,hpoints);

xi = linspace(basis.knot{1}(1),basis.knot{1}(end),shp(2));

% fig = PlotRefinement(basis,hspace,X0,Y0,Z0);

figure
plot(xi,hfuncs0')

%% Refinement loop

RefFunc = [1 7 1];

i = 0;
for iFunc = RefFunc
    
    i = i + 1;
    
    disp('____________')
    disp(['Refinement ',num2str(i)])
    
    hspace  = RefineHspace(basis, hspace, iFunc);
    hpoints = RefineHpoints(hspace, hpoints);
    
    disp(['Hierarchcal level ', num2str(hspace.level)])
    disp(['Number of functions ', num2str(hspace.nfunc)])
    
    [X,Y,~,hfuncs] = ConstructGeometry(basis,hspace,hpoints);
    
    xi = linspace(basis.knot{1}(1),basis.knot{1}(end),shp(2));
    
%     fig = PlotRefinement(basis,hspace,X,Y,Z,hpoints);
%     title(['Refinement ',num2str(i)])
%     view(-45,30)

    figure
    plot(xi,hfuncs')
    
    figure
    plot(X,Y)

end

