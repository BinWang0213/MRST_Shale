close all;clear;

%Langmuir adsorption
load('Case1_Langmuir_10000mdft_PROS.mat');
load('Case1_Langmuir_GridEDFM_PROS_EDFM.mat');
load('Case1_Langmuir_10000mdft_PROS_EDFM.mat');

Data_Ref1 = csvread('CMG_Desorption10000mdft.csv',6,1);

%% Gas flow rate Plot
%MRST-Shale
t0=time_list;
y0=Case1_Fcd10000mdft_GasPro;
y1=Case1_Grid_EDFM_GasProEDFM; 
y2=Case1_Fcd10000mdft_GasProEDFM; 

%CMG
t4=Data_Ref1(:,1); y4=Data_Ref1(:,3).*ft^3; 
logidx4=round(exp(linspace(log(1),log(numel(t4)),65)));


figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;


l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b--', 'LineWidth', 3);
l3=plot(t0, y2,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
plot(t4(logidx4), y4(logidx4),'ko', 'LineWidth', 1,'MarkerSize',9);
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
legend('LGR', ...
       'EDFM',...
       'LGR+EDFM',...
       'Commerical simulator',...
       'FontSize',20)
   
%% Cumulative Gas flow rate Plot
Cum0=cumtrapz(t0,y0); 
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);

Cum4=cumtrapz(t4,y4);

Diff_EDFM=(Cum1(end)-Cum0(end))/Cum0(end);
Diff_LGREDFM=(Cum2(end)-Cum0(end))/Cum0(end);


logidx4=round(exp(linspace(log(1),log(numel(t4)),150)));

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0/365, Cum0/1e6,'r-', 'LineWidth', 3);
l2=plot(t0/365, Cum1/1e6,'b--', 'LineWidth', 3);
l4=plot(t0/365, Cum2/1e6,'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
plot(t4/365, Cum4/1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
hold off;

grid on;
box on;

set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production [10^6 m^3]'))
xlim([0 30]);
legend('LGR', ...
       'EDFM +3.04%',...
       'LGR+EDFM -0.27%',...
       'Commerical simulator',...
       'FontSize',20)