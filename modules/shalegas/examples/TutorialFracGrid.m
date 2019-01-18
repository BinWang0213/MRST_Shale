%% Tutorial of grid generation

close all; clear;
mrstModule add shalegas

%% Problem 1 Single fracture with Log refinement
% Create explicit fracture grid with Log LGR
physdim = [1990*ft 1990*ft];
DZ=150*ft;

%Define Fracture geometry
NumFracs=1;
Frac_Spacing=1990*ft;
Frac_halfLength=350*ft;
Frac_height=DZ; %Thickness of reservoir
Frac_StartXY=[physdim(1)/2 physdim(2)/2-Frac_halfLength];
NumHydraulicFracs=1;
NumNaturalFracs=0;

%Create Planar Hydraulic Fractures
[fl,xy_wells]=createMultiStageFracs(NumFracs,Frac_Spacing,...
    Frac_halfLength,Frac_StartXY);


%Example 1. Coarse 2D example
%%{
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
    'NX_FracRefine',5,... % NX for each half fracture spacing
    'NY_FracRefine',2,...% NY for half fracture length
    'NY_OutRefine',1,... % NY for region out of fracture
    'FracCellSize',10*ft,... % DX for fracture cell, fracture aperature 
    'FracCellSize_Y',35*ft); %DY for fracture cell
%%}
G = computeGeometry(G);
[NX,NY]=deal(G.cartDims(1),G.cartDims(2));

%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',600,'Title','Example 1');
%Plot Grid
plotGrid(G,'FaceAlpha',0.5,'EdgeAlpha',0.8), view(2), axis equal tight;


%Example 2. Fine 3D example
physdim = [1990*ft 1990*ft 100*ft]; % 3D Grid requires 3D physical dimension
%%{
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
    'NX_FracRefine',50,... % NX for each half fracture spacing
    'NY_FracRefine',5,...% NY for half fracture length
    'NY_OutRefine',3,... % NY for region out of fractured zone 
    'FracCellSize',10*ft,... % DX for fracture cell, fracture aperature 
    'FracCellSize_Y',35*ft); %DY for fracture cell
%%}

%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',600,'Title','Example 2');
%Plot Grid
plotGrid(G,'FaceAlpha',0.5,'EdgeAlpha',0.8), view(3), axis equal tight;


%% Problem 2 Multiple Fracture with Log refinement
physdim = [200 140 10];

%Define Fracture geometry
NumFracs=5;
Frac_Spacing=25;
Frac_halfLength=30;
Frac_StartXY=[50 40];
NumHydraulicFracs=5;
NumNaturalFracs=0;

%Planar Hydraulic Fractures
[fl,xy_wells]=createMultiStageFracs(NumFracs,Frac_Spacing,...
    Frac_halfLength,Frac_StartXY);
Total_length=calcFracsLength(fl)


%Option2. Local grid refinement for vertical hydraulic fractures
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
   'NX_FracRefine',10,...  % NX for each half fracture spacing
   'NX_OutRefine',2,...    % NX for region out of fractured zone 
   'NY_FracRefine',5,...% NY for half fracture length
   'NY_OutRefine',2,... % NY for region out of fractured zone 
   'FracCellSize',1e-3,...
   'FracCellSize_Y',1.5763);

G = computeGeometry(G);
[NX,NY]=deal(G.cartDims(1),G.cartDims(2));

%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',600,'Title','Example 3');
%Plot Grid
plotGrid(G,'FaceAlpha',0.5,'EdgeAlpha',0.8),view(3), axis equal tight;


