%% Benchmark problem from http://dx.doi.org/10.1016/j.jngse.2015.08.013

close all; clear;
%PathConfigure;
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-props ad-core % AD framework
mrstModule add ad-blackoil      % Three phase simulator
mrstmodule add OpenShale

%% Define geometric quantitites
physdim = [200 140];
DZ=10*meter;

%Define fracture position
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

%Non-planar 1 Hydraulic Fractures
%{
rng(123456790); %Fix the random pattern
[fl,xy_wells]=createMultiStageNonPlanarFracs(NumFracs,Frac_Spacing,...
    Frac_halfLength,Frac_StartXY,'Severity',0.3,'ShapeFactor',1);
Total_length=calcFracsLength(fl)
%}

%Non-planar 2 Hydraulic Fractures
%{
[fl1,len1]=createNonPlanarFrac(xy_wells(1,:),0.25,Frac_halfLength*1.1,'left');
[fl2,len2]=createNonPlanarFrac(xy_wells(2,:),0.5,Frac_halfLength*0.6,'left');
[fl4,len4]=createNonPlanarFrac(xy_wells(4,:),0.5,Frac_halfLength*0.6,'right');
[fl5,len5]=createNonPlanarFrac(xy_wells(5,:),0.25,Frac_halfLength*1.1,'right');
[fl3,len3]=createPlanarFrac(xy_wells(3,:),(Total_length-len1-len2-len4-len5)/2);
fl=[fl1;fl2;fl3;fl4;fl5];
Total_length=calcFracsLength(fl) %Check total length
%}

%% Define Grid

%Option1. Unifrom mesh for irregular fractures
%celldim = [187 183];
%G = cartGrid(celldim, physdim);

%Option2. Local grid refinement for vertical hydraulic fractures
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
   'NX_FracRefine',20,'NX_OutRefine',8,...
   'FracCellSize',0.001,'FracCellSize_Y',1.5763);

G = computeGeometry(G);
G.NumNFs=NumNaturalFracs;
G.NumHFs=NumHydraulicFracs;


%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',600,'Title','EDFM Grid');
%Plot Grid
plotGrid(G,'FaceAlpha',1,'EdgeAlpha',0.1), box on,view(2), axis equal tight;

%% Process fracture lines
[G,fracture] = processFracture2D(G,fl);
fracture.aperture = 1e-3*meter; % Fracture aperture
fracture.height=DZ; %Thickness of reservoir
G.fractureHeight=fracture.height;
clf; plotFractureLines(G,fracture);axis equal tight;


%% Compute CI and construct fracture grid
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
%cell_size = physdim(1)/G.cartDims(1); min_size = cell_size/2;  % minimum and average cell size.
cell_size = mean(physdim)/mean(G.cartDims); min_size = cell_size/2;  % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); axis equal tight;box on;

%% Define rock properties
p0=160*barsa;
G.rock.perm = ones(G.cells.num,1) * 100*nano*darcy;
G.rock.poro = ones(G.cells.num,1) * 0.1;
K_frac = 1;% Darcy in EDFM
%Fracture conductivity
FractureConductivity=K_frac*1000*convertTo(fracture.aperture,ft);
Fcd=K_frac*fracture.aperture/(G.rock.perm(1)*Frac_halfLength);
G = makeRockFrac(G, K_frac,'porosity', 1-1e-8);

[G,T] = defineNNCandTrans(G,F,fracture); %Define nnc trans

%% Black-oil shale gas fluid properties
fluid=setShaleGasFluid_Case2(G,G.rock,p0);

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
damage_perm=0.0037*G.rock.perm(1);
S=pi/2*damage_extent/Frac_halfLength*(G.rock.perm(1)/damage_perm-1);
W = addWellEDFM([], G, G.Matrix.rock, cellInx, ...
        'Type', 'bhp', 'Val', 40*barsa, ...
        'Radius', 0.1,'Comp_i',1,'Skin',S);

%time step has to be setup with wells
numSteps = 1;                  % number of time-steps
totTime  = 10000*day;
dt       = totTime / numSteps;  % constant time step
%dt_list=repmat(dt,numSteps,1);
rampSteps=27;
dt_list = rampupTimesteps(totTime, dt,rampSteps); %simple rampup steps
time_list=cumsum(convertTo(dt_list,day));

schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
p0=repmat(p0, [G.cells.num, 1]);
%p0(G.Matrix.cells.num+1:end)=40*barsa;
state  = initResSol(G, p0);%0-single phase model

%% Run simulations
[ws_e, states_e, report_e] = simulateScheduleAD(state, model, schedule);

if isfield(fluid,'mG_ad')
    data_file='Jiang2015_FullMechs.csv';
else
    data_file='Jiang2015_NoMechs.csv';
end

PlotEDFMGasRate(time_list,ws_e, ...
    'Formation_thickness',DZ,... %meters
    'Reference_data',data_file,...
    'YUnit', meter^3/day,...
    'XUnit', day,...
    'Xlim',[1e-4 1e4],...
    'Ylim',[1e2 1e5],...
    'CumPlot',1,...
    'LogLog',1);
PlotEDFMPresSurf(fl,G,states_e,numel(time_list),'ColorLim',[4 16])
