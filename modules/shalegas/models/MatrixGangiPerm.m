function [k_gangi]=MatrixGangiPerm(FracCellMask,parameters)
%{
Set up matrix perm reduction due to micro-fracture clousure

Reference
---------
Gangi,1978 10.1016/0148-9062(78)90957-9
Shale data1 SPE-18267-PA
Shale data2 10.1016/j.coal.2015.12.014 table 4

Arguments
---------
Pc    -- Overburden confining pressure, Pa
m     -- Gangi surface roughness coefficient, dimensionless
P1    -- matrix effective stress when micro-fracture fully closed, Pa
alpha -- biot conefficient
p     -- pore pressure,Pa

Return
---------
[k_gangi] permeability corrfection term due to fracture clousure

Author: Bin Wang(binwang.0213@gmail.com)
Date: Dec.2018
%}


m=0.5; 
Pc=38*mega*Pascal; %1.1 psi/feet GEOMECHANICAL STUDIES OF THE BARNETT SHALE, TEXAS, USA
P1=180*mega*Pascal;
alpha=0.5;

k_gangi_func=@(p) (1-((Pc-alpha.*p)./P1).^m).^3;

%Natural fracture closure only occurs on the matrix
k_gangi = @(p) ~FracCellMask.*k_gangi_func(p)+FracCellMask;
%k_gangi = @(p) FracCellMask.*k_gangi_func(p)+~FracCellMask;

%{
%Verify with Gangi (1978) Figure 8
p    = linspace(0*Pascal, 20*1e3*Pascal, 2000);
m=0.2;
k_gangi_test=@(p) (1-((p)./max(p)).^m).^3;
figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
plot(p./max(p), k_gangi_test(p),'bo-', 'LineWidth', 2)
set(gca,'FontSize',25);
xlim([0,0.04])
xlabel('Normalized Pressure')
ylabel('Gangi F_{frac}')
title('Gangi Perm Correction Factor','FontSize',20)


%p    = linspace(40*barsa, 140*barsa, 50);
p    = linspace(0, 300*barsa, 50);
figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
%plot((Pc-alpha.*p)./P1, k_gangi_func(p),'bo-', 'LineWidth', 2)
plot(convertTo(p,barsa), k_gangi_func(p),'bo-', 'LineWidth', 2)
set(gca,'FontSize',25);
grid on;
ylim([0.15 0.3])
xlabel('Pressure [Bar]')
ylabel('F_{app}')
title('Gangi Perm Correction Factor','FontSize',20)

%}

end