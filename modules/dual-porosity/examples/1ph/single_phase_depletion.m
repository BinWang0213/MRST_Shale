%% This is a single-phase dual-porosity example
% We compare the results of a two simulations: one with two layers of cells 
% (one layer representing the fracture and the other the matrix) and a
% second simulation with a single layer of cells but with dual-porosity
% transfer. The only transfer mechanism is depletion. Synce there is just 
% one layer in the matrix, the results should match perfectly. 

%% Read the problem from a deckfile
mrstModule add ad-props ad-core ad-blackoil dual-porosity

%% Create grid 
% Two layer grid
x_size_2l = 500;
y_size_2l = 500;
z_size_2l = 10;

Nx_2l = 50;
Ny_2l = 2;
Nz_2l = 1;

G_2l = cartGrid([Nx_2l Ny_2l Nz_2l],[x_size_2l y_size_2l z_size_2l]);
G_2l = computeGeometry(G_2l);

% One layer grid (Dual-Porosity)
x_size_1l = 500;
y_size_1l = 250;
z_size_1l = 10;

Nx_1l = 50;
Ny_1l = 1;
Nz_1l = 1;

dx_1l = x_size_1l/Nx_1l;
dy_1l = y_size_1l/Ny_1l;
dz_1l = z_size_1l/Nz_1l;

G_1l = cartGrid([Nx_1l Ny_1l Nz_1l],[x_size_1l y_size_1l z_size_1l]);
G_1l = computeGeometry(G_1l);

%% Create rock
% Base porosity and permeability
km = 0.1*milli*darcy;
kf = 1000*km;
phi = 0.2;

% Two layer grid
rock = makeRock(G_2l, km, phi);
bottom_layer = find(G_2l.cells.centroids(:,2)<=y_size_2l/2);
rock.perm(bottom_layer,:) = kf;

% One layer grid (Dual-Porosity)
rock_fracture = makeRock(G_1l, kf, phi);
rock_matrix = makeRock(G_1l, km, phi);

%% Define pressures (drawdown)
pres = 8000*psia;
pwf = 1000*psia;

%% Create fluid
% Same fluid properties for both continua
fluid = initSimpleADIFluid('mu', [0.1 1], 'rho', [1 1 1], 'n', [2 2 2]);

b = 1/1.2; %reciprocal of formation vol factor.
c = 1e-07/barsa;
fluid.bW = @(p) b*exp((p - pres)*c);    

fluid_fracture = fluid;
fluid_matrix = fluid;

%% Define Wells 
W_2l = verticalWell([],G_2l,rock,1,1,'type','bhp','Val',pwf,'Comp_i',1);
W_1l = verticalWell([],G_1l,rock_fracture,1,1,'type','bhp','Val',pwf,'Comp_i',[1,0]);

dt = repmat(1*day, 20, 1);

schedule_2l = simpleSchedule(dt, 'W', W_2l);
schedule_11 = simpleSchedule(dt, 'W', W_1l);
% Set up reservoir
gravity off

%% Shape factor
As = dx_1l*dz_1l;
Ls = dy_1l/2;
Vb = dx_1l*dy_1l*dz_1l;

shape_factor = (As/Ls)*(1/Vb);

%% Create the models
model_2l = WaterModel(G_2l, rock, fluid);
model_1l = TwoPhaseOilWaterDPModel(G_1l, rock_fracture, fluid_fracture, rock_matrix, fluid_matrix, []);

% Here is the important bit for a dual-porosity simulation: you have to set
% the transfer_model_object to one of the transfer functions in the folder 
% transfer_models. Here we use the the same transfer function that is defined
% in SLB's Eclipse simulator. Note that we also define a shape factor value 
% by setting the shape_factor_value field of the object shape_factor_object
% inside the transfer model.
model_1l.transfer_model_object = EclipseTransferFunction();
model_1l.transfer_model_object.shape_factor_object.shape_factor_value = shape_factor;

%% Initializing state
state_2l.pressure = ones(G_2l.cells.num,1)*pres;
state_2l.s = repmat([1 0],G_2l.cells.num,1);
state_2l.wellSols= initWellSolAD(W_2l, model_2l, state_2l);

state_1l.pressure = ones(G_1l.cells.num,1)*pres;
state_1l.pom = state_1l.pressure;
state_1l.s = repmat([1 0],G_1l.cells.num,1);
state_1l.swm = ones(G_1l.cells.num,1);
state_1l.wellSols= initWellSolAD(W_1l, model_1l, state_1l);

%% Simulate the models
[~, states_2l] = simulateScheduleAD(state_2l, model_2l, schedule_2l);
[~, states_1l] = simulateScheduleAD(state_1l, model_1l, schedule_11);


%% Plotting Information
pressure_cells = find(G_2l.cells.centroids(:,2)>=y_size_2l/2);
ti = 0;

fig1 = figure('Position',[100,100,900,400]);
for i = 1:length(states_1l)
    ti = ti + dt(i);

    figure(fig1)
    subplot(1,2,1)
    title(['Pressure @ T = ' num2str(ti/day) ' days'])
    p = plotCellData(G_2l,states_2l{i}.pressure/psia);
    plotWell(G_2l,W_2l)

    axis equal tight off
    view(-20, 50)
    colorbar
    caxis([pwf/psia pres/psia])

    figure(fig1)
    subplot(1,2,2)
    plot(G_1l.cells.centroids(:,1),states_1l{i}.pom(:,1)/psia,'r+');
    hold on
    plot(G_2l.cells.centroids(pressure_cells,1),states_2l{i}.pressure(pressure_cells,1)/psia,'LineWidth',1.5,'LineStyle','-','Color','b');
    hold off
    grid on
    ylim([pwf/psia pres/psia])
    xlabel('x [m]')
    ylabel('P [psia]')
    title(['Pressure in the matrix @ T = ' num2str(ti/day) ' days'])
    legend({'Dual Porosity','2 Layer Model'},'Location','Southeast')

    drawnow;
end














