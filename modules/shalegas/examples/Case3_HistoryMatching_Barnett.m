%% Case3. Benchmark problem from https://doi.org/10.1016/j.fuel.2016.03.055

close all; clear;
%PathConfigure;
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-props ad-core % AD framework
mrstModule add ad-blackoil      % Three phase simulator
mrstModule add shalegas

%% Define geometric quantitites
physdim = [1200 300];
DZ=90*meter;

%Define fracture position
NumFracs=28;
Frac_Spacing=30.5;
Frac_halfLength=47.2;
Frac_StartXY=[188.333 102.333];
NumHydraulicFracs=5;
NumNaturalFracs=0;

%Planar Hydraulic Fractures
[fl,xy_wells]=createMultiStageFracs(NumFracs,Frac_Spacing,...
    Frac_halfLength,Frac_StartXY);
Total_length=calcFracsLength(fl)

%% Define Grid

%Option1. Unifrom mesh for irregular fractures
%celldim = [3001 31];
%G = cartGrid(celldim, physdim);

%Option2. Local grid refinement for vertical hydraulic fractures
%On
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
   'NX_FracRefine',4,'NX_OutRefine',5,'FracCellSize',0.1);

G = computeGeometry(G);
G.NumNFs=NumNaturalFracs;
G.NumHFs=NumHydraulicFracs;


%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',300,'Title','Fracture Map');
%Plot Grid
plotGrid(G,'FaceAlpha',1,'EdgeAlpha',0.3),view(2), axis tight;

%% Process fracture lines
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 3e-3*meter; % Fracture aperture
fracture.height=DZ; %Thickness of reservoir
G.fractureHeight=fracture.height;
clf; plotFractureLines(G,fracture);axis equal tight;

%% Compute CI and construct fracture grid
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
cell_size = mean(physdim)/mean(G.cartDims); min_size = cell_size/2;  % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); axis equal tight;box on;

%% Define rock properties
p0=203.4*barsa;
G.rock.perm = ones(G.cells.num,1) * 150*nano*darcy;
G.rock.poro = ones(G.cells.num,1) * 0.03;
K_frac = 0.1;% Darcy in EDFM
%Fracture conductivity
FractureConductivity=K_frac*1000*convertTo(fracture.aperture,ft);
Fcd=K_frac*fracture.aperture/(G.rock.perm(1)*Frac_halfLength);
G = makeRockFrac(G, K_frac,'porosity', 1-1e-8);

[G,T] = defineNNCandTrans(G,F,fracture); %Define nnc trans

%% Black-oil shale gas fluid properties
fluid=setShaleGasFluid_Case3(G,G.rock,p0);

%% Define shale gas flow model
model = ShaleGasModel(G,[],fluid);

N = getNeighbourship(G, 'topological', true);
intx = all(N ~= 0, 2);

% Send in internal transmissibility and neighborship to the operator setup
model.operators = setupOperatorsTPFA(G, G.rock, 'trans', T(intx), 'neighbors', N(intx, :));

%% Assume constant BHP horizontal well
[cellInx]=markFracWellCell(xy_wells,G);

%Correction Skin - using face damage skin equation
%http://www.fekete.com/san/webhelp/welltest/webhelp/Content/HTML_Files/Reference_Materials/Skin.htm
damage_extent=0.1*Frac_halfLength; % meter
damage_perm=0.008*G.rock.perm(1);
S=pi/2*damage_extent/Frac_halfLength*(G.rock.perm(1)/damage_perm-1);
W = addWellEDFM([], G, G.Matrix.rock, cellInx, ...
        'Type', 'bhp', 'Val', 34.5*barsa, ...
        'Radius', 0.1,'Comp_i',1,'Skin',S);

%time step has to be setup with wells
numSteps = 1;                  % number of time-steps
%totTime  = 1600*day;
totTime  = 30*year;
dt       = totTime / numSteps;  % constant time step
%dt_list=repmat(dt,numSteps,1);
rampSteps=25;
dt_list = rampupTimesteps(totTime, dt,rampSteps); %simple rampup steps
time_list=cumsum(convertTo(dt_list,day));

schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
p0=repmat(p0, [G.cells.num, 1]);
state  = initResSol(G, p0);%0-single phase model

%% Run simulations
[ws_e, states_e, report_e] = simulateScheduleAD(state, model, schedule);

data_file='Barnett_PRO_Field.csv';
PlotEDFMGasRate(time_list,ws_e, ...
    'Formation_thickness',DZ,... %meters
    'Reference_data',data_file,...
    'YUnit', meter^3/day,...
    'XUnit', day,...
    'Xlim',[1e-1 1e4],...
    'Ylim',[1e4 1e6],...
    'LogLog',0);
PlotEDFMPresSurf(fl,G,states_e,numel(time_list),'ColorLim',[])
%PlotEDFMGasRateTransient(time_list,ws_e,'Formation_thickness',fracture.height);