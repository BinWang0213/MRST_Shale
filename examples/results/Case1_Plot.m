close all;clear;
load('Case1_Langmuir_PROS_EDFM.mat');
load('Case1_base_PROS_EDFM.mat');
%load('Case1_Langmuir_PROS.mat');
%load('Case1_base_PROS.mat');
Data_Ref = csvread('CMG_PRO_base.csv',1);
Data_Ref2 = csvread('CMG_PRO_Langmuir.csv',1);

%% Gas flow rate Plot
%MRST-Shale
t0=time_list;
y0=Case1_Langmuir_GasProEDFM_MRST_Shale;
y1=Case1_base_GasProEDFM_MRST_Shale;
%y0=Case1_Langmuir_GasPro_MRST_Shale;
%y1=Case1_base_GasPro_MRST_Shale;
%CMG GEM
t2=Data_Ref2(:,1); y2=Data_Ref2(:,2).*ft^3; 
t3=Data_Ref(:,1); y3=Data_Ref(:,2).*ft^3;

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0, y0,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0, y1,'b--', 'LineWidth', 3);
logidx=round(exp(linspace(log(1),log(numel(t3)),40)));
l3=plot(t2(logidx), y2(logidx),'ko', 'LineWidth', 1,'MarkerSize',9);
l4=plot(t3(logidx), y3(logidx),'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;

grid on;
set(gca,'FontSize',20);
xlabel('Time [Days]')
ylabel('Gas Flow Rate [m^3/day]')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
xlim([1e-4 1e4]);
%ylim([1e4 1e7]);
legend('MRST-Shale with adsorpton', ...
       'MRST-Shale without adsorpton',...
       'Commerical simulator')

%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t2,y2);
Cum3=cumtrapz(t3,y3);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0/365, Cum0./1e6,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0/365, Cum1./1e6,'b--', 'LineWidth', 3);
l3=plot(t2(1:5:end)/365, Cum2(1:5:end)./1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
l4=plot(t3(1:5:end)/365, Cum3(1:5:end)./1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;

set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production\n[10^6 m^3]'))
xlim([0 30]);
legend('MRST-Shale with adsorpton', ...
       'MRST-Shale without adsorpton',...
       'Commerical simulator')

%% Pressure over line plot
P_dist_CMG = csvread('CMG_Pressure_overLine.csv',2);
P_dist_MRST_Shale = csvread('Case1_MRST_Shale_Pressure_overline.csv',2);

d0=P_dist_CMG(:,1)*ft;p0=P_dist_CMG(:,2)*psia;
d1=P_dist_MRST_Shale(:,1)*ft;p1=P_dist_MRST_Shale(:,2)*psia;

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(d1, convertTo(p1,barsa),'r-', 'LineWidth', 3)
hold on;
l2=plot(d0, convertTo(p0,barsa),'b--', 'LineWidth', 3)
grid on;
hold off;

set(gca,'FontSize',20);
xlabel('Distance [m]')
ylabel('Pressure [Bar]')
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
xlim([0 850]);
ylim([0 400]);
legend('MRST-Shale', ...
       'Commerical simulator')

