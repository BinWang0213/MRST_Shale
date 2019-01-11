close all;clear;
PathConfigure;

%Langmuir adsorption
load('Case1_Langmuir_GridLGR_NFs_PROS.mat');
load('Case1_Langmuir_GridEDFM_NFs_PROS_EDFM.mat');
load('Case1_Langmuir_GridLGREDFM_NFs_PROS_EDFM.mat');

load('Case1_Langmuir_GridLGREDFM_PROS_EDFM.mat');

%% Gas flow rate Plot
%MRST-Shale
t0=time_list;
y0=Case1_GridEDFM_NFs_GasProEDFM; 
y1=Case1_GridLGR_NFs_GasPro;
y2=Case1_GridLGREDFM_NFs_GasProEDFM;

%CMG
t4=t0; y4=Case1_Grid_LGREDFM_GasProEDFM; 
logidx4=round(exp(linspace(log(1),log(numel(t4)),60)));


figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;


l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
l3=plot(t0, y2,':', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
%l4=plot(t4(logidx4), y4(logidx4),'ko', 'LineWidth', 1,'MarkerSize',9);
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
%ylim([1e3 1e8]);
legend('LGR with NFs', ...
       'EDFM with NFs',...
       'LGR+EDFM with NFs',...
       'LGR without NFs',...
       'FontSize',20)
   
%% Cumulative Gas flow rate Plot
Cum0=cumtrapz(t0,y0); 
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);

Cum4=cumtrapz(t4,y4);

Diff_EDFM=(Cum1(end)-Cum0(end))/Cum0(end);
Diff_LGREDFM=(Cum2(end)-Cum0(end))/Cum0(end);


logidx4=round(exp(linspace(log(1),log(numel(t4)),120)));

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0/365, Cum0/1e6,'r-', 'LineWidth', 3);
l2=plot(t0/365, Cum1/1e6,'b--', 'LineWidth', 3);
l4=plot(t0/365, Cum2/1e6,':', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
%l4=plot(t4/365, Cum4/1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;

grid on;
box on;

set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production [10^6 m^3]'))
xlim([0 30]);
legend('LGR with NFs', ...
       'EDFM with NFs +1.58%',...
       'LGR+EDFM with NFs +1.49%',...
       'LGR without NFs',...
       'FontSize',20)
