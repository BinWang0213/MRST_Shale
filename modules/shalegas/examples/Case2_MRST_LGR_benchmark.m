%% Benchmark problem from http://dx.doi.org/10.1016/j.jngse.2015.08.013

close all; clear;
%PathConfigure;
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-props ad-core % AD framework
mrstModule add ad-blackoil      % Three phase simulator
mrstModule add shalegas

%% Define geometric quantitites
physdim = [200 140 10];

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

%% Define Grid

%Option1. Unifrom mesh for irregular fractures
%celldim = [187 183];
%G = cartGrid(celldim, physdim);

%Option2. Local grid refinement for vertical hydraulic fractures
G = ExplicitFracGrid(physdim,...
    NumFracs,Frac_Spacing,Frac_halfLength,Frac_StartXY,...
   'NX_FracRefine',20,'NX_OutRefine',8,...
   'FracCellSize',1e-3,'FracCellSize_Y',1.5763);

G = computeGeometry(G);
[NX,NY]=deal(G.cartDims(1),G.cartDims(2));
G.NumNFs=NumNaturalFracs;
G.NumHFs=NumHydraulicFracs;


%Plot Fracture Map
plotFracGeo(physdim,fl,xy_wells,'FigSize',600);
%Plot Grid
plotGrid(G,'FaceAlpha',0.5,'EdgeAlpha',0.1), box on,view(2), axis equal tight;

%% Define rock properties
p0=160*barsa;
[poro_rock,poro_frac]=deal(0.1,1.0);
[perm_rock,perm_frac]=deal(100*nano*darcy,1*darcy);
gravity reset off;

perm=repmat(perm_rock,NX,NY);
perm(G.FracCell.I,G.FracCell.J)=perm_frac;
poro=repmat(poro_rock,NX,NY);
poro(G.FracCell.I,G.FracCell.J)=poro_frac;

rock = makeRock(G, perm(:), poro(:));

%% Black-oil shale gas fluid properties
fluid=setShaleGasFluid_Case2(G,rock,p0);

%% Define shale gas flow model
model = WaterModelG(G,rock,fluid);

%% Assume constant BHP horizontal well
IJ_wells = markCellbyXY(xy_wells,G);
cellInx = sub2ind(G.cartDims, IJ_wells(:,1), IJ_wells(:,2));
W = addWell([], G, rock, cellInx,'Dir', 'x','Radius', 0.1, ...
        'Type', 'bhp', 'Val', 40*barsa,'Comp_i',1);

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


%% Plot well and permeability
% Since the well is inside the reservoir, we remove a section around the
% well so that we can see the well path
%{
clf;
show = true([G.cells.num, 1]);
show(cellInx) = false;% Hide well cell
plotCellData (G , convertTo(rock.perm,milli*darcy),show, 'EdgeColor', 'k');
%plotCellData(G, convertTo(p_init, barsa), show, 'EdgeColor', 'k');
colorbar ('horiz'); view (2); axis equal tight;
%}

%% Run simulations
[ws_e, states_e, report_e] = simulateScheduleAD(state, model, schedule);

if isfield(fluid,'mG_ad')
    data_file='Jiang2015_FullMechs.csv';
else
    data_file='Jiang2015_NoMechs.csv';
end

PlotEDFMGasRate(time_list,ws_e, ...
    'Reference_data',data_file,...
    'YUnit', meter^3/day,...
    'XUnit', day,...
    'Xlim',[1e-4 1e4],...
    'Ylim',[1e2 1e5],...
    'LogLog',1);
PlotEDFMPresSurf(fl,G,states_e,numel(time_list),'ColorLim',[4 16])
%PlotEDFMGasRateTransient(time_list,ws_e,'Formation_thickness',fracture.height);