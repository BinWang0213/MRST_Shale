%% Case 1 Benchmark with CMG GEM
% In this example, we will demonstrate how one can easily extend the
% compressible single-phase pressure solver to include the effect of
% pressure-dependent viscosity using either arithmetic averaging of the
% viscosity or harmonic averaging of the fluid mobility.
% DataFolder: '\examples\Benchmark_CMG\'
clear;close all;clc;
% mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core ad-props mrst-gui compositional

% ---------------------------------------------------------------------------
% Basic parameters setup
%---------------------------------------------------------------------------

%% Define geometric quantitites
% Create explicit fracture grid with Log LGR
physdim = [1990*ft 1990*ft 150*ft];

%Define fracture geometry
NumFracs=1;
Frac_Spacing=500*ft;
Frac_halfLength=350*ft;
Frac_height=150*ft; %Thickness of reservoir
Frac_StartXY=[physdim(1)/2 physdim(2)/2-Frac_halfLength];

[fl,xy_wells]=createMultiStageFracs(NumFracs,Frac_Spacing,...
    Frac_halfLength,Frac_StartXY);

%% Define geometric quantitites
% Create explicit fracture grid with Log LGR
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
    'NX_FracRefine',100,...
    'NX_OutRefine',15,...
    'NY_FracRefine',15,...
    'NY_OutRefine',10,...
    'FracCellSize',0.01*ft,...
    'FracCellSize_Y',0.01*ft);
G = computeGeometry(G);

[NX,NY]=deal(G.cartDims(1),G.cartDims(2));

%Plot Grid
% plotGrid(G,'FaceAlpha',1.0,'EdgeAlpha',0.4), view(2), axis equal tight;

%% Define rock properties
[poro_rock,poro_frac]=deal(0.07,1.0);
[perm_rock,perm_frac]=deal(0.0005*milli*darcy,0.5*darcy);


perm=repmat(perm_rock,NX,NY);
perm(G.FracCell.I,G.FracCell.J)=perm_frac;
poro=repmat(poro_rock,NX,NY);
poro(G.FracCell.I,G.FracCell.J)=poro_frac;

rock = makeRock(G, perm(:), poro(:));

%% Compositional shale gas fluid properties
% Name of problem and pressure range
casename = 'onlymethane';
pwf = 500*psia;
p0 = 5000*psia;

%% Assume constant BHP vertical well
IJ_wells = markCellbyXY(xy_wells,G);
cellInx = sub2ind(G.cartDims, IJ_wells(:,1), IJ_wells(:,2));

W = [];
% Producer
W = addWell(W, G, rock, cellInx, 'Dir', 'x','Radius', 0.25*ft, ...
    'comp_i', [1], 'Val', pwf, 'sign', -1, 'Type', 'bhp');    


% ---------------------------------------------------------------------------
% Compositional model
%---------------------------------------------------------------------------
%% Set up model and initial state
nkr = 2;
[fluid, info] = getCompositionalFluidCase(casename);
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);

gravity reset off
model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', false);
% model = OverallCompositionCompositionalModel(G, rock, flowfluid, fluid, 'water', false);

ncomp = fluid.getNumberOfComponents();
s0 = [1];
TinK = 327.594;
state0 = initCompositionalState(G, p0, TinK, s0, info.initial, model.EOSModel);

for i = 1:numel(W)
    W(i).components = info.initial;
end
%% Set up schedule and simulate the problem
%time step has to be setup with wells
M = csvread('CMG_timestep2.csv',1);
dt_list=M(:,1)*day;
time_list=cumsum(convertTo(dt_list,day));

schedule = simpleSchedule(dt_list, 'W', W);

%% Run simulations
[ws_comp, states_comp, report_comp] = simulateScheduleAD(state0, model, schedule);

%---------------------------------------------------------------------------
% Black oil model
%---------------------------------------------------------------------------
%% Black-oil shale gas fluid properties
[fluid]=setShaleGasFluid_Case1comp(G,rock);

%% Define shale gas flow model
model = ShaleGasModel(G,rock,fluid);
schedule = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state0  = initResSol(G, p0, 0);%0-single phase model

%% Run simulations
[ws_BO, states_BO, report_BO] = simulateScheduleAD(state0, model, schedule);



%% Comparsion figures
names = {'Compositional', 'BlackOil'};
ws = {ws_comp, ws_BO};
shortname = {'comp', 'BO'};
plotWellSols(ws, cumsum(schedule.step.val), 'datasetnames', names)


figure
plotToolbar(G, states_comp);
axis equal tight off
daspect([1 1 0.2])
view(0, 90);
plotWell(G, W);
colormap jet(25)
title(names{1});
colorbar('horiz')

figure; plotToolbar(G, states_BO);
axis equal tight off
daspect([1 1 0.2])
view(0, 90);
plotWell(G, W);
colormap jet(25)
title(names{2});
colorbar('horiz')


PV = sum(G.cells.volumes .* rock.poro)/(ft^3)  %41582 Mrcf
% Bgi = (14.6959/520)*z*TinK/p0;  %rcf/scf
% STOIIP = sum(G.cells.volumes .* rock.poro)/Bgi; %in scm
