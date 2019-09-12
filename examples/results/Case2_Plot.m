close all;clear;
load('Case2_FullMechs_PROS_EDFM.mat');
load('Case2_Langmuir_PROS_EDFM.mat');
load('Case2_Slippage_PROS_EDFM.mat');
load('Case2_NoMechs_PROS_EDFM.mat');
Data_Ref = csvread('Jiang2015_FullMechs.csv',1);
Data_Ref2 = csvread('Jiang2015_NoMechs.csv',1);

%% Gas flow rate Plot
%MRST-Shale
t0=time_list;
y0=Case2_FullMechs_GasProEDFM_MRST_Shale;
y1=Case2_Langmuir_GasProEDFM_MRST_Shale;
y2=Case2_Slippage_GasProEDFM_MRST_Shale;
y3=Case2_NoMechs_GasProEDFM_MRST_Shale;

t4=Data_Ref2(:,1); y4=Data_Ref2(:,2); 
t5=Data_Ref(:,1); y5=Data_Ref(:,2);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
l3=plot(t0, y2,':', 'LineWidth', 3,'color',[43/255,140/255,190/255]);
l4=plot(t0, y3,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
l5=plot(t4, y4,'ko', 'LineWidth', 1,'MarkerSize',9);
l6=plot(t5, y5,'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;

grid on;
set(gca,'FontSize',20);
xlabel('Time [Days]')
ylabel('Gas Flow Rate [m^3/day]')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
box on;
xlim([1e-4 1e4]);
ylim([1e1 1e5]);
legend('Full Mechanism', ...
       'Adsorption only',...
       'Slippage & diffusion only',...
       'No Mechanism',...
       'Jiang and Younis, 2015')

%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);
Cum3=cumtrapz(t0,y3);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0/365, Cum0./1e6,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0/365, Cum1./1e6,'b--', 'LineWidth', 3);
l2=plot(t0/365, Cum2./1e6,':', 'LineWidth', 3,'color',[43/255,140/255,190/255]);
l2=plot(t0/365, Cum3./1e6,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
%l3=plot(t2(1:5:end)/365, Cum2(1:5:end)./1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
%l4=plot(t3(1:5:end)/365, Cum3(1:5:end)./1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;
grid on;


set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production\n[10^6 m^3]'))
xlim([0 25]);
ylim([0 5]);
legend('Full Mechanism', ...
       'Adsorption only',...
       'Slippage & diffusion only',...
       'No Mechanism')
   
%% Irregular fracture shape
clear;
load('Case2_FullMechs_Planar_PROS_EDFM.mat');
load('Case2_FullMechs_NonPlanar1_PROS_EDFM.mat');
load('Case2_FullMechs_NonPlanar2_PROS_EDFM.mat');

t0=time_list;
y0=Case2_Planar_GasProEDFM_MRST_Shale;
y1=Case2_NonPlanar1_GasProEDFM_MRST_Shale;
y2=Case2_NonPlanar2_GasProEDFM_MRST_Shale;

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;
l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
l3=plot(t0, y2,':', 'LineWidth', 3,'color',[43/255,140/255,190/255]);

hold off;

grid on;
set(gca,'FontSize',20);
xlabel('Time [Days]')
ylabel('Gas Flow Rate [m^3/day]')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
box on;
xlim([1e-4 1e4]);
ylim([1e1 1e5]);
legend('Planar', ...
       'Non-planar1',...
       'Non-planar2')
   
%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);

Diff_1=(Cum1(end)-Cum0(end))/Cum0(end);
Diff_2=(Cum2(end)-Cum0(end))/Cum0(end);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0/365, Cum0./1e6,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0/365, Cum1./1e6,'b--', 'LineWidth', 3);
l2=plot(t0/365, Cum2./1e6,':', 'LineWidth', 3,'color',[43/255,140/255,190/255]);
%l3=plot(t2(1:5:end)/365, Cum2(1:5:end)./1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
%l4=plot(t3(1:5:end)/365, Cum3(1:5:end)./1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;
grid on;


set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production\n[10^6 m^3]'))
xlim([0 25]);
ylim([0 5]);
legend('Planar', ...
       'Non-planar1',...
       'Non-planar2')