close all;clear;
load('Case3_HistoryMatching_PROS_EDFM');
Data_Ref = csvread('Barnett_PRO_Field.csv',1);

%% History matching
%MRST-Shale
t0=time_list;
y0=Case3_HistoryMatching_GasProEDFM;

%Field Data (Cao, 2016)
t4=Data_Ref(:,1); y4=Data_Ref(:,2);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l2=plot(t0, y0,'r-', 'LineWidth', 3);
l3=plot(t4, y4,'ko', 'LineWidth', 1,'MarkerSize',9);
l3.Color(4) = 0.2;

hold off;

grid on;
box on;
set(gca,'FontSize',20);
xlabel('Time [Days]')
ylabel('Gas Flow Rate [m^3/day]')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
xlim([1e-0 1e4]);
ylim([1e4 1e6]);
legend('MRST-Shale', ...
       'Field Data')

%% Production forcast

load('Case3_Forecast_PROS_EDFM.mat');
%MRST-Shale
t0=time_list;
y0=Case3_Forecast_GasProEDFM;

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');


l1=plot(t0./365, y0,'r-', 'LineWidth', 3)
hold on;
l3=plot(t4/365, y4,'ko', 'LineWidth', 1,'MarkerSize',9);
l3.Color(4) = 0.2;

grid on;

hold off;

set(gca,'FontSize',20);
xlabel('Time [Year]')
ylabel('Gas Flow Rate [m^3/day]')
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
xlim([0 30]);
ylim([0 2e5]);
legend('MRST-Shale', ...
       'Field Data')