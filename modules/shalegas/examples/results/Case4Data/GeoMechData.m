close all;clear;
RefData = csvread('ExperimenalFcd.csv',3);



%% propped fractures
p_propped=RefData(1:11,9)*psia; 
FcdN_propped_stiff=RefData(1:11,10);
FcdN_propped_med=RefData(1:11,12);
FcdN_propped_soft=RefData(1:11,14);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

S_close = linspace(300*psia, 10000*psia, 20);

l1=plot(convertTo(S_close,barsa), 10.^(-0.00011.*convertTo(S_close,psia)-0.0971),'r-', 'LineWidth', 3);
l2=plot(convertTo(S_close,barsa), 10.^(-0.00036.*convertTo(S_close,psia)+0.2396),'b--', 'LineWidth', 3);
l3=plot(convertTo(S_close,barsa), 10.^(-0.00064.*convertTo(S_close,psia)+0.4585),'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
plot(convertTo(p_propped,barsa), FcdN_propped_stiff,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(convertTo(p_propped,barsa), FcdN_propped_med,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(convertTo(p_propped,barsa), FcdN_propped_soft,'ko', 'LineWidth', 1,'MarkerSize',9);

hold off;

grid on;
box on;
set(gca,'FontSize',20);
xlabel('Closure pressure [barsa]')
ylabel('Normalized Fracture Conductivity [-]')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
%xlim([1e-2 1e4]);
%ylim([1e3 1e6]);
legend('Stiff Shale', ...
       'Median Shale', ...
       'Soft Shale',...
       'Alramahi and Sundberg (2012)')

%% unpropped fractures
p_unpropped_stiff=RefData(1:8,1)*psia; FcdN_unpropped_stiff=RefData(1:8,2);
p_unpropped_med=RefData(1:5,3)*psia; FcdN_unpropped_med=RefData(1:5,4);
p_unpropped_soft=RefData(1:5,5)*psia; FcdN_unpropped_soft=RefData(1:5,6);


figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
hold on;

S_close = linspace(300*psia, 10000*psia, 20);

l1=plot(convertTo(S_close,barsa), 10.^(-0.793.*log(convertTo(S_close,psia))+4.5618),'r-', 'LineWidth', 3);
l2=plot(convertTo(S_close,barsa), 10.^(-0.89.*log(convertTo(S_close,psia))+5.0725),'b--', 'LineWidth', 3);
l3=plot(convertTo(S_close,barsa), 10.^(-1.041.*log(convertTo(S_close,psia))+6.0216),'-.', 'LineWidth', 3,'color',[44/255, 162/255, 95/255]);
plot(convertTo(p_unpropped_stiff,barsa), FcdN_unpropped_stiff,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(convertTo(p_unpropped_med,barsa), FcdN_unpropped_med,'ko', 'LineWidth', 1,'MarkerSize',9);
plot(convertTo(p_unpropped_soft,barsa), FcdN_unpropped_soft,'ko', 'LineWidth', 1,'MarkerSize',9);

hold off;

grid on;
box on;
set(gca,'FontSize',20);
xlabel('Closure pressure [barsa]')
ylabel('Normalized Fracture Conductivity [-]')
set(gca, 'YScale', 'log')

set(gca, 'MinorGridColor', 'w')
%xlim([1e-2 1e4]);
%ylim([1e3 1e6]);
legend('Stiff Shale', ...
       'Median Shale', ...
       'Soft Shale',...
       'Wu et al (2018)')