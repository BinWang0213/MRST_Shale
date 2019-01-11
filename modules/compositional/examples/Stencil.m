%% Example demonstrating a three dimensional, six component problem
% We set up a simple grid and rock structure. The numbers can be adjusted
% to get a bigger/smaller problem. The default is a small problem with
% 20x20x2 grid blocks.
close all

mrstModule add ad-core ad-props mrst-gui compositional
% Dimensions

% Name of problem and pressure range
casename = 'onlymethane';
pwf = 500*psia;
p_init = 5000*psia;

% Set up grid and rock
x_max = 120;
x_min = 0;
max_nLogGrids  = 100;
deltaR         = 0.05;  % Frac width will be 2*deltaR    0.001meters = 1mm 
factor         = 1.3;

pt = zeros(1,10);
pt(1) = deltaR;
fprintf('Pt[%d]= %24.16f \n',1,pt(1));
for ilogSpacnNum = 2:max_nLogGrids 
   deltaR = deltaR * factor;
   pt(ilogSpacnNum) = pt(ilogSpacnNum-1) + deltaR;
   if(pt(ilogSpacnNum) >= ((x_max-x_min)/2.0))
      fprintf('num of log grids is: %d \n', ilogSpacnNum);
      break
   end
   fprintf('Pt[%d]= %24.16f \n',ilogSpacnNum,pt(ilogSpacnNum));
end

ilogSpacnNum = ilogSpacnNum - 1;
pt = pt(1:ilogSpacnNum);
myX = 60 + [-fliplr(pt) ,pt];
%increase from 5 to 20
temp = 320+cumsum(fliplr([6.3956  8.3643   10.9236   14.2507   18.5759]));

myY = [0:20:320 temp 385:5:500];
G = tensorGrid(myX, myY,[0 50]);
G.cells.num
plotGrid(G,'FaceColor' ,[ .7 .7 1]); 

G = computeGeometry(G);
% rock = makeRock(G, 100*nano*darcy, 0.04);
rock = makeRock(G, 1.0e-19, 0.04);
%set fracture properties
fracIndx_Y = 28:G.cartDims(2);
numFracCells = length(fracIndx_Y);
fracIndx_X = ones(numFracCells,1)*ceil(G.cartDims(1)/2);
fraccells = sub2ind(G.cartDims, fracIndx_X, fracIndx_Y');
rock.poro(fraccells) = 0.33*3/100;
rock.perm(fraccells) = 5e-11; %50*darcy;

%% Set up quarter five spot well pattern
% We place vertical wells in opposing corners of the reservoir. The
% injector is rate controlled and the producer is bottom hole pressure
% controlled.
W = [];
% Producer
ind = [ceil(G.cartDims(1)/2) G.cartDims(2)];
prodcell = sub2ind(G.cartDims, ind(1), ind(2));
W = verticalWell(W, G, rock, prodcell, [],'Radius', 0.1, ...
    'comp_i', [1], 'Name', 'Prod', 'Val', pwf, 'sign', -1, 'Type', 'bhp');
%% Set up model and initial state
% We set up a problem with quadratic relative permeabilities. The fluid
% model is retrieved from "High Order Upwind Schemes for Two-Phase,
% Multicomponent Flow" (SPE 79691) by B. T. Mallison et al.
%
% The model consists of six components. Several of the components are not
% distinct molecules, but rather lumped groups of hydrocarbons with similar
% molecular weight. The reservoir contains all these components initially,
% which is then displaced by the injection of pure CO2.

nkr = 2;
[fluid, info] = getCompositionalFluidCase(casename);
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);

gravity reset on
model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', false);
% model = OverallCompositionCompositionalModel(G, rock, flowfluid, fluid, 'water', false);

ncomp = fluid.getNumberOfComponents();
s0 = [1];
state0 = initCompositionalState(G, p_init, info.temp, s0, info.initial, model.EOSModel);

for i = 1:numel(W)
    W(i).components = info.initial;
end
%% Set up schedule and simulate the problem
% We simulate two years of production with a geometric rampup in the
% timesteps.
time = 30*year;
dt = rampupTimesteps(time, 2*year, 30);
schedule = simpleSchedule(dt, 'W', W);

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%% Plot all the results
lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, -200, 350, 200]);
nm = ceil(ncomp/2);
v = [-30, 60];
for step = 1:numel(states)
    figure(h); clf
    state = states{step};
    for i = 1:ncomp
        subplot(nm, 3, i);
        plotCellData(G, state.components(:, i), 'EdgeColor', 'none');
        view(v);
        title(fluid.names{i})
        caxis([0, 1])
    end
    subplot(nm, 3, ncomp + 1);
    plotCellData(G, state.pressure, 'EdgeColor', 'none');
    view(v);
    title('Pressure')
    
%     subplot(nm, 3, ncomp + 2);
%     plotCellData(G, state.s(:, 1), 'EdgeColor', 'none');
%     view(v);
%     title('sO')
    
%     subplot(nm, 3, ncomp + 3);
%     plotCellData(G, state.s(:, 2), 'EdgeColor', 'none');
%     view(v);
%     title('sG')
%     drawnow
end
%% Plot the results in the interactive viewer
figure(1); clf;
plotToolbar(G, states)
view(v);
axis tight


figure(2);
plotWellSols(ws,cumsum(schedule.step.val)/86400)
tinDays = cumsum(schedule.step.val)/86400;

%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
