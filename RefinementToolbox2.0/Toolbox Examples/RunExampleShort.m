clear all, close all

%% Load NURBS data

load('ExamSurf')
nurbs = ExamSurf.srf;

%% Settings

truncated = 1;          
n         = 200;
dim       = 3;
reffunc   = [1 9 12 15 18 21 24];

%% Initialization

basis   = MakeBasis(nurbs, n);
hspace  = MakeHspace(nurbs, truncated);
hpoints = MakeHpoints(nurbs, dim);


for i = reffunc
%% Refinement 
tstart = tic;
hspace  = RefineHspace(basis, hspace, i);
hpoints = RefineHpoints(hspace, hpoints);
toc(tstart);
%% Geometry construction

end 

[X,Y,Z] = ConstructGeometry(basis, hspace, hpoints);

%% Plotting

fig = PlotRefinement(basis, hspace, X, Y, Z, hpoints);