close all;clear;
load('Case4_Planar_PROS_EDFM.mat');
load('Case4_Nonplanar_PROS_EDFM.mat');
load('Case4_Nonplanar_NFs_PROS_EDFM.mat');
load('Case4_Nonplanar_NFs_Geomechs_PROS_EDFM.mat');


%% Gas rate

%MRST-Shale
t0=time_list;
y0=Case4_Planar_GasProEDFM;
y1=Case4_Nonplanar_GasProEDFM;
y2=Case4_Nonplanar_NFs_GasProEDFM;

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
l3=plot(t0, y2,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);

hold off;

grid on;
box on;
set(gca,'FontSize',20);
xlabel('Time [Days]')
ylabel('Gas Flow Rate [m^3/day]')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
xlim([1e-2 1e4]);
ylim([1e3 1e6]);
legend('Planar', ...
       'Non-planar', ...
       'Non-planar + NFs')


%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);

%relative error 
Diff_1=(Cum1(end)-Cum0(end))/Cum0(end);
Diff_2=(Cum2(end)-Cum0(end))/Cum0(end);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0/365, Cum0./1e6,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0/365, Cum1./1e6,'b--', 'LineWidth', 2.9);
l3=plot(t0/365, Cum2./1e6,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
hold off;
grid on;


set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production\n[10^6 m^3]'))
xlim([0 30]);
ylim([0 150]);
legend('Planar (reference)', ...
       'Non-planar -5.69%', ...
       'Non-planar + NFs,  +8.63%')

%% Geomechanics history matching
Data_Ref = csvread('Barnett_PRO_Field.csv',1);

t0=time_list;
y0=Case4_Planar_GasProEDFM;
y1=Case4_Nonplanar_NFs_Geomechs_GasProEDFM;

%Field Data (Cao, 2016)
t4=Data_Ref(:,1); y4=Data_Ref(:,2);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
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
xlim([1e-2 1e4]);
ylim([1e3 1e6]);
legend('Planar', ...
       'Non-planar+NFs+Geomechs', ...
       'Field Data')

%% Production forcast

%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);

%relative error 
Diff_3=(Cum1(end)-Cum0(end))/Cum0(end);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0/365, Cum0./1e6,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0/365, Cum1./1e6,'b--', 'LineWidth', 2.9);
hold off;
grid on;


set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production [10^6 m^3]'))
xlim([0 30]);
ylim([0 150]);
legend('Planar', ...
       'Non-planar+NFs+Geomechs -0.97%')
