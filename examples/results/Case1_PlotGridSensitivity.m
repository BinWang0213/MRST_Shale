close all;clear;
Data_Ref1 = csvread('Case1_LGR50.csv',1);
Data_Ref2 = csvread('Case1_LGR150.csv',1);
Data_Ref3 = csvread('Case1_LGR250.csv',1);
Data_Ref4 = csvread('Case1_LGR300.csv',1);

%% Gas flow rate Plot
%MRST-Shale
t0=Data_Ref1(:,1);%day
y0=Data_Ref1(:,2);%m3/day
y1=Data_Ref2(:,2);
y2=Data_Ref3(:,2);
y3=Data_Ref4(:,2);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0, y0,'r-', 'LineWidth', 3);
l2=plot(t0, y1,'b-', 'LineWidth', 3);
l4=plot(t0, y2,'-', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
l3=plot(t0, y3,'k:', 'LineWidth', 3);
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
%ylim([1e4 1e7]);
legend('NX-FracRefine=50', ...
       'NX-FracRefine=150',...
       'NX-FracRefine=250',...
       'NX-FracRefine=300')

%% Cumulative production rate
Cum0=cumtrapz(t0,y0);
Cum1=cumtrapz(t0,y1);
Cum2=cumtrapz(t0,y2);
Cum3=cumtrapz(t0,y3);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

l1=plot(t0/365, Cum0/1e6,'r-', 'LineWidth', 3);
l2=plot(t0/365, Cum1/1e6,'b-', 'LineWidth', 3);
l4=plot(t0/365, Cum2/1e6,'-', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
l3=plot(t0/365, Cum3/1e6,'k:', 'LineWidth', 3);
hold off;

grid on;
box on;
set(gca,'FontSize',20);
xlabel('Time [Years]')
ylabel(sprintf('Cumulative Gas Production\n[10^6 m^3]'))
xlim([0 30]);
legend('NX-FracRefine=50', ...
       'NX-FracRefine=150',...
       'NX-FracRefine=250',...
       'NX-FracRefine=300')


