%% quarterSpotDP - This example outlines the basic usage of the 
% two-phase dual-porosity model. We inject the first phase on 
% one corner and produce at the opposite one.
clc;
close all;

%% Include the essential modules. dual-porosity model has the functionality
% we need.
mrstModule clear
mrstModule add ad-props ad-core ad-blackoil blackoil-sequential dual-porosity

%% Create the grid: a simple cartesian grid, 100m x 100m and 50 blocks on x
% and y directions.

% Two layer grid
x_size = 100;
y_size = 100;
z_size = 1;

Nx = 20;
Ny = 20;
Nz = 1;

G = cartGrid([Nx Ny Nz],[x_size y_size z_size]);
G = computeGeometry(G);

%% We need two rock structures for the dual-porosity model: one for the fracture
% system and one for the non-fractured rock matrix. Fractures usually will 
% have high permeability and low pore volume, while matrix usually will 
% have a higher porosity and decreased permeability. Therefore,
% interconnected fractures generally serves as a "highway" for fluid flow, 
% while the fluid transfer with the matrix will happen at a larger
% timescale
kf = 10000*milli*darcy;
phif = 0.01;

km = 1*milli*darcy;
phim = 0.1;

% Rock Structures
rock_fracture = makeRock(G, kf, phif);
rock_matrix = makeRock(G, km, phim);

%% Pressures of injector and producer wells
pres = 0*psia;
pprod = pres - 1000*psia;
pinj = pres + 1000*psia;

%% We also need two fluid structures for the dual-porosity model. Fractures are
% usually set as having zero capillary pressure, but we don't define any 
fluid_matrix = initSimpleADIFluid('mu', [1 2], 'rho', [1 1 1], 'n', [2 2 1]);
fluid_fracture = fluid_matrix;

b = 1/1.2; %reciprocal of formation vol factor.
c = 1e-07/barsa;
fluid_fracture.bW = @(p) b*exp((p - pres)*c);    
fluid_fracture.bO = @(p) b*exp((p - pres)*c);    

%% We set pc = 0 in the fractures and a linear 
% capillary pressure in the matrix. When injecting water in a fractured 
% reservoir, the main mechanim of oil recovery is spontaneous imbibition
% of the water in the fractures into the matrix.
Pcscale = 50*kilo*Pascal;
fluid_fracture.pcOW=@(swm)0;
fluid_matrix.pcOW=@(swm)-Pcscale*swm + Pcscale;

%% Add the wells
W = addWell([],G,rock_fracture,Nx*Ny,'type','bhp','Val',pprod,'Comp_i',[1,1]);
W = addWell(W,G,rock_fracture,1*1,'type','bhp','Val',pinj,'Comp_i',[1,0]);

%% Create the model. TwoPhaseOilWaterDPModel is a modified version of 
% TwoPhaseOilWaterModel that adds dual porosity behaviour
gravity off

%% Model definition
model = TwoPhaseOilWaterDPModel(G, rock_fracture, fluid_fracture,...
                               rock_matrix, fluid_matrix, []);
                               
% The shape factor and transfer function
fracture_spacing = repmat([1,1,1],G.cells.num,1);
shape_factor_name = 'KazemiShapeFactor';
model.transfer_model_object = EclipseTransferFunction(shape_factor_name,fracture_spacing);

%% Initialize the field with a constant pressure and fully saturated by oil
state0.pressure = ones(G.cells.num,1)*pres;
state0.s = repmat([0 1],G.cells.num,1);
state0.swm = zeros(G.cells.num,1);
state0.pom = ones(G.cells.num,1)*pres;

%% Initialize the well solution
state0.wellSol= initWellSolAD(W, model, state0);
state = state0;

solver = NonLinearSolver();

%% Handles to pictures that are going to be plotted at each timestep
fig1 = figure(1);

%% Source
src_val = 0.0001*sum(poreVolume(G,rock_fracture))/day;
src = addSource([], 200, src_val, 'sat', [0,1]);

%% Simulate the models
dT = 0.1*day;
n = 20;
for i = 1:n
    %% Advancing fields
    state = solver.solveTimestep(state, dT, model, 'W', W, 'src', src);
    disp(i)

    %% Plotting fields
    figure(fig1)
    subplot(2,2,1)
    title(sprintf('Sw in the fractures'))
    p = plotCellData(G,state.s(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)

    figure(fig1)
    subplot(2,2,2)
    title(sprintf('Sw in the matrix'))
    p = plotCellData(G,state.swm(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)

    figure(fig1)
    subplot(2,2,3)
    title(sprintf('Water transfer rate'))
    p = plotCellData(G,state.Twm(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)

    figure(fig1)
    subplot(2,2,4)
    title(sprintf('Oil transfer rate'))
    p = plotCellData(G,state.Tom(:,1));
    p.EdgeAlpha = 0;
    axis equal tight off
    view(-20, 50)

    drawnow;
end
