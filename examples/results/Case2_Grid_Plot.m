close all;clear;
load('Case2_Grid161x39_LGR_PROS_EDFM.mat');
load('Case2_Grid161x39_Uniform_PROS_EDFM.mat');
load('Case2_Grid1001x39_Uniform_PROS_EDFM.mat');
Data_Ref2 = csvread('Jiang2015_NoMechs.csv',1);

%% Gas flow rate Plot
%MRST-Shale
t0=time_list;
y0=Case2_Grid161x39_GasProEDFM;
y1=Case2_UniformGrid161x39_GasProEDFM;
y2=Case2_UniformGrid1001x39_GasProEDFM;

t4=Data_Ref2(:,1); y4=Data_Ref2(:,2); 

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
l4=plot(t0, y2,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
l5=plot(t4, y4,'ko', 'LineWidth', 1,'MarkerSize',9);
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
legend('LGR 161x39',...
       'Uniform 161x39',...
       'Uniform 1001x39',...
       'Jiang and Younis, 2015')

%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');

l1=plot(t0/365, Cum0./1e6,'r-', 'LineWidth', 3);
hold on;
l2=plot(t0/365, Cum1./1e6,'b--', 'LineWidth', 3);
l3=plot(t0/365, Cum2./1e6,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
hold off;
grid on;


set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production\n[10^6 m^3]'))
xlim([0 25]);
ylim([0 2]);
legend('LGR 161x39',...
       'Uniform 161x39',...
       'Uniform 1001x39')
   
