close all;clear;
%PathConfigure;
%{ 
%Non-desorption
load('Case1_base_10000mdft_PROS.mat');
load('Case1_base_10000mdft_PROS_EDFM.mat');
load('Case1_base_50mdft_PROS.mat');
load('Case1_base_50mdft_PROS_EDFM.mat');
load('Case1_base_5mdft_PROS.mat');
load('Case1_base_5mdft_PROS_EDFM.mat');

Data_Ref1 = csvread('CMG_base10000mdft.csv',6,1);
Data_Ref2 = csvread('CMG_base50mdft.csv',6,1);
Data_Ref3 = csvread('CMG_base5mdft.csv',6,1);
%}

%Langmuir adsorption
load('Case1_Langmuir_10000mdft_PROS.mat');
load('Case1_Langmuir_10000mdft_PROS_EDFM.mat');
load('Case1_Langmuir_50mdft_PROS.mat');
load('Case1_Langmuir_50mdft_PROS_EDFM.mat');
load('Case1_Langmuir_5mdft_PROS.mat');
load('Case1_Langmuir_5mdft_PROS_EDFM.mat');

Data_Ref1 = csvread('CMG_Desorption10000mdft.csv',6,1);
Data_Ref2 = csvread('CMG_Desorption50mdft.csv',6,1);
Data_Ref3 = csvread('CMG_Desorption5mdft.csv',6,1);

%% Gas flow rate Plot
%MRST-Shale
t0=time_list;
y0=Case1_Fcd10000mdft_GasPro; y0e=Case1_Fcd10000mdft_GasProEDFM;
y1=Case1_Fcd50mdft_GasPro; y1e=Case1_Fcd50mdft_GasProEDFM;
y2=Case1_Fcd5mdft_GasPro; y2e=Case1_Fcd5mdft_GasProEDFM;

%CMG
t3=Data_Ref1(:,1); y3=Data_Ref1(:,3).*ft^3; 
t4=Data_Ref2(:,1); y4=Data_Ref2(:,3).*ft^3; 
t5=Data_Ref3(:,1); y5=Data_Ref3(:,3).*ft^3; 
logidx3=round(exp(linspace(log(1),log(numel(t3)),55)));
logidx4=round(exp(linspace(log(1),log(numel(t4)),33)));
logidx5=round(exp(linspace(log(1),log(numel(t5)),39)));

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;


l1=plot(t0, y0e,'r-', 'LineWidth', 3);
l2=plot(t0, y1e,'b-', 'LineWidth', 3);
l4=plot(t0, y2e,'-', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
plot(t3(logidx3), y3(logidx3),'ko', 'LineWidth', 1,'MarkerSize',9);
plot(t0, y0,'r--', 'LineWidth', 3);
plot(t4(logidx4), y4(logidx4),'ko', 'LineWidth', 1,'MarkerSize',9);
plot(t0, y1,'b--', 'LineWidth', 3);
plot(t5(logidx5), y5(logidx5),'ko', 'LineWidth', 1,'MarkerSize',9);
plot(t0, y2,'--', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
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
ylim([1e3 1e8]);
legend('Fcd=10000 md{\cdot}ft', ...
       'Fcd=50 md{\cdot}ft',...
       'Fcd=5 md{\cdot}ft',...
       'Commerical simulator',...
       'FontSize',20)
   
%% Cumulative Gas flow rate Plot
Cum0=cumtrapz(t0,y0); Cum0e=cumtrapz(t0,y0e); 
Cum1=cumtrapz(t0,y1);Cum1e=cumtrapz(t0,y1e);
Cum2=cumtrapz(t0,y2);Cum2e=cumtrapz(t0,y2e);

Cum3=cumtrapz(t3,y3);
Cum4=cumtrapz(t4,y4);
Cum5=cumtrapz(t5,y5);

Diff_10000mdft=(Cum0e(end)-Cum0(end))/Cum0(end);
Diff_50mdft=(Cum1e(end)-Cum1(end))/Cum1(end);
Diff_5mdft=(Cum2e(end)-Cum2(end))/Cum2(end);


logidx3=round(exp(linspace(log(1),log(numel(t3)),60)));
logidx4=round(exp(linspace(log(1),log(numel(t4)),60)));
logidx5=round(exp(linspace(log(1),log(numel(t5)),70)));

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0/365, Cum0e/1e6,'r-', 'LineWidth', 3);
l2=plot(t0/365, Cum1e/1e6,'b-', 'LineWidth', 3);
l4=plot(t0/365, Cum2e/1e6,'-', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
plot(t3/365, Cum3/1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(t0/365, Cum0/1e6,'r--', 'LineWidth', 3);
plot(t4/365, Cum4/1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(t0/365, Cum1/1e6,'b--', 'LineWidth', 3);
plot(t5/365, Cum5/1e6,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(t0/365, Cum2/1e6,'--', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
hold off;

grid on;
box on;

set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production [10^6 m^3]'))
xlim([0 30]);
legend('Fcd=10000 md{\cdot}ft', ...
       'Fcd=50 md{\cdot}ft',...
       'Fcd=5 md{\cdot}ft',...
       'Commerical simulator',...
       'FontSize',20)
