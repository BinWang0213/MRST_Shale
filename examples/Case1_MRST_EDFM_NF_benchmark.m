%% Pressure-dependent viscosity
% In this example, we will demonstrate how one can easily extend the
% compressible single-phase pressure solver to include the effect of
% pressure-dependent viscosity using either arithmetic averaging of the
% viscosity or harmonic averaging of the fluid mobility.
% DataFolder: '\examples\Benchmark_CMG\'

close all; clear;
%PathConfigure;
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-props ad-core % AD framework
mrstModule add ad-blackoil      % Three phase simulator
mrstModule add OpenShale

%% Define geometric quantitites
% Create explicit fracture grid with Log LGR
physdim = [1990*ft 1990*ft];
DZ=150*ft;


%Define fracture geometry
NumFracs=1;
Frac_Spacing=500*ft;
Frac_halfLength=350*ft;
Frac_height=150*ft; %Thickness of reservoir
Frac_StartXY=[physdim(1)/2 physdim(2)/2-Frac_halfLength];
NumHydraulicFracs=1;
NumNaturalFracs=0;

[fl,xy_wells]=createMultiStageFracs(NumFracs,Frac_Spacing,...
    Frac_halfLength,Frac_StartXY);


%% Define Grid

%Option1. Uniform grid
%celldim = [201 159];
%G = cartGrid(celldim, physdim);

%Option2. Local grid refinement for vertical hydraulic fractures

%{ 
%EDFM
[G,fl] = ExplicitFracGridNF(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
    'NX_FracRefine',100,'NX_OutRefine',15,...
    'NY_OutRefine',10,...
    'FracCellSize',0.01*ft,...
    'FracCellSize_Y',6*ft,...
    'EDFM_Grid',1,...
    'NumNFs',3,'NY_Refine',21,'NY_LogRefine',false,...
    'NF_Spacing',2*Frac_halfLength/3,...
    'NF_Length',249.7*ft,...
    'NF_RepeatPatternSpace',0*ft,...
    'NF_StartXY',[Frac_StartXY(1)-249.7*ft Frac_StartXY(2)]);
%}

%%{ 
%EDFM + LGR
[G,fl] = ExplicitFracGridNF(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
    'NX_FracRefine',100,'NX_OutRefine',15,...
    'NY_OutRefine',10,...
    'FracCellSize',0.01*ft,...
    'FracCellSize_Y',0.01*ft,...
    'EDFM_Grid',1,...
    'NumNFs',3,'NY_FracRefine',25,'NY_LogRefine',true,...
    'NF_Spacing',2*Frac_halfLength/3,...
    'NF_Length',249.7*ft,...
    'NF_RepeatPatternSpace',0*ft,...
    'NF_StartXY',[Frac_StartXY(1)-249.7*ft Frac_StartXY(2)]);
%%}
NumHydraulicFracs=1;
NumNaturalFracs=6;

G = computeGeometry(G);
[NX,NY]=deal(G.cartDims(1),G.cartDims(2));
G.NumNFs=NumNaturalFracs;
G.NumHFs=NumHydraulicFracs;

%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',600,'Title','LGR+EDFM Grid');
%Plot Grid
plotGrid(G,'FaceAlpha',0.3,'EdgeAlpha',0.3), view(2), axis equal tight;

%% Process fracture lines
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 0.01*ft; % Fracture aperture
fracture.height = Frac_height; %Thickness of reservoir
G.fractureHeight=fracture.height;
clf; plotFractureLines(G,fracture);axis equal tight;

%% Compute CI and construct fracture grid
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
cell_size = max(physdim)/max(G.cartDims)*1.1; min_size = cell_size/2.1;  % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); axis equal tight;box on;

%% Define rock properties
p0=5000*psia;
G.rock.perm = ones(G.cells.num,1) * 500*nano*darcy;
G.rock.poro = ones(G.cells.num,1) * 0.07;
K_frac = 1000;% [Darcy] unit in EDFM
%Fracture conductivity
FractureConductivity=K_frac*1000*convertTo(fracture.aperture,ft); %md-ft
Fcd=K_frac*fracture.aperture/(G.rock.perm(1)*Frac_halfLength);
G = makeRockFrac(G, K_frac,'porosity', 1-1e-8);

%Set natural fractures permeability
K_NF=0.5; %darcy
G = updateRockFrac(G, K_NF,'NaturalFracs');

[G,T] = defineNNCandTrans(G,F,fracture); %Define nnc trans

%% Black-oil shale gas fluid properties
[fluid]=setShaleGasFluid_Case1(G,G.rock);

%% Define shale gas flow model
model = ShaleGasModel(G,[],fluid);

N = getNeighbourship(G, 'topological', true);
intx = all(N ~= 0, 2);

% Send in internal transmissibility and neighborship to the operator setup
model.operators = setupOperatorsTPFA(G, G.rock, 'trans', T(intx), 'neighbors', N(intx, :));

%% Assume constant BHP horizontal well
[cellInx]=markFracWellCell(xy_wells,G);
W = addWellEDFM([], G, G.Matrix.rock, cellInx, 'Dir', 'x',...
        'Type', 'bhp', 'Val', 500*psia, ...
        'Radius', 0.25*ft,'Comp_i',1);

%time step has to be setup with wells
M = csvread('CMG_timestep2.csv',1);%Benchmark_CMG\CMG_timestep.csv
dt_list=M(:,1)*day;
time_list=cumsum(convertTo(dt_list,day));

schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
p0=repmat(p0, [G.cells.num, 1]);
state  = initResSol(G, p0, 0);%0-single phase model

%% Run simulations
[ws_e, states_e, report_e] = simulateScheduleAD(state, model, schedule);

PlotEDFMGasRate(time_list,ws_e, ...
    'Formation_thickness',DZ,... %meters
    'YUnit', meter^3/day,...
    'XUnit', day,...
    'Xlim',[1e-4 1e4],...
    'CumPlot',1,...
    'LogLog',1);
PlotEDFMPresSurf(fl,G,states_e,numel(time_list),'ColorLim',[3.5 35])
