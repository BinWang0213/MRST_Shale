%% EXPERIMENTAL: Set up and run a five-spot problem with water-alternating-gas 
% (WAG) drive. One approach to hydrocarbon recovery is to inject gas or water. In the
% water-alternating-gas approach, wells change between injection of gas and
% water to improve the sweep efficiency.

clear;
clc;
close all;

% We begin by loading the required modules
mrstModule clear
mrstModule add ad-core ad-blackoil ad-props dual-porosity

[ state, model, schedule ] = three_phase_water_alternating_gas_setup();

%% Simulate the schedule
[~, states] = simulateScheduleAD(state, model, schedule);

%% Plotting results
fig = figure(1);
for i = 1:numel(states)
    
    figure(fig)
    subplot(2,3,1)
    title('Sw - matrix')
    p = plotCellData(model.G,states{i}.swm);
    p.EdgeAlpha = 0;
    grid on
    
    figure(fig)
    subplot(2,3,2)
    title('So - matrix')
    p = plotCellData(model.G,1-states{i}.swm-states{i}.sgm);
    p.EdgeAlpha = 0;
    grid on
    
    figure(fig)
    subplot(2,3,3)
    title('Sg - matrix')
    p = plotCellData(model.G,states{i}.sgm);
    p.EdgeAlpha = 0;
    grid on
    
    figure(fig)
    subplot(2,3,4)
    title('Sw - fractures')
    p = plotCellData(model.G,states{i}.s(:,1));
    p.EdgeAlpha = 0;
    grid on
    
    figure(fig)
    subplot(2,3,5)
    title('So - fractures')
    p = plotCellData(model.G,1-states{i}.s(:,1)-states{i}.s(:,3));
    p.EdgeAlpha = 0;
    grid on
    
    figure(fig)
    subplot(2,3,6)
    title('Sg - fractures')
    p = plotCellData(model.G,states{i}.s(:,3));
    p.EdgeAlpha = 0;
    grid on
    
    drawnow
   
end

