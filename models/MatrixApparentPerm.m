function [k_app]=MatrixApparentPerm(method,rock,mu,FracCellMask,parameters)
%{
Set up apparent permeability term

viscous flow (Kn ? 0.001), slip flow (0.001 ? Kn ? 0.1),
transition flow (0.1 ? Kn ? 10) and Knudsen flow (Kn ? 10)

Arguments
---------
method -- 'Langmuir' 'BET'
Kn     -- Knudsen number,dimensionless
alpha  -- rarefaction parameter, dimensionless


Return
---------
[m_ad] [func of accumation term considering adsorption]

Author: Bin Wang(binwang.0213@gmail.com)
Date: Nov.2018
%}


if(strcmp(method,'Civan'))
    MW=parameters.MW; %kg/mol
    R=8.314;   %J/(K.mol)
    T=parameters.T_Res;  %K 343.15
    
    Kn      = @(p) mu(p)./2.8284./p.*sqrt(3.141592653*R*T/2/MW.*rock.poro./rock.perm);
    %This is used for MRST-AD, MRST-AD doesn't support sin,cos,atan,....
    %Taylor expansion can be used instead
    alpha   = @(p) 128/15/pi^2.*atan(4.*Kn(p).^0.4); %Added h=atan in ADI.m
    %alpha   = @(p) 1.2;
    %This is exact for plot
    %alpha   = @(p) 128/15/pi^2.*atan(4.*Kn(p).^0.4);
    k_app_func   = @(p) (1+alpha(p).*Kn(p)).*(1+4.*Kn(p)./(1+Kn(p)));

end

%Gas slippage and diffusion only occures on tight matrix
k_app = @(p) ~FracCellMask.*k_app_func(p)+FracCellMask;

%{

%Real Knudsen number
p    = linspace(1*barsa, 300*barsa, 100);
left_color = [0 0 1];
right_color = [1 0 0];
figure('rend','painters','pos',[10 10 800 600], ...
       'defaultAxesColorOrder',[left_color; right_color]);
Kn_sample=Kn(p);

yyaxis left
plot(convertTo(p, barsa), Kn_sample(1,:),'b-', 'LineWidth', 2)
set(gca, 'YScale', 'log')
ylabel('Kn')
ylim([1e-1 1e1]);

yyaxis right
k_app_sample=k_app(p);
plot(convertTo(p, barsa), k_app_sample(1,:),'r--', 'LineWidth', 2)
ylabel('Fapp')
ylim([1e0 1e2]);

set(gca, 'YScale', 'log')
set(gca,'FontSize',25);
xlabel('Pressure [Bar]')
grid on;
title('Perm Correction factor(Florence et al, 2007)','FontSize',20)

%Matched well with 10.1007/s11242-011-9842-6
figure('rend','painters','pos',[10 10 1200 600]);
Kn    = linspace(1e-3, 1e2, 2000);
alpha   = @(Kn) 128/15/pi^2.*atan(4.*Kn.^0.4);
k_app   = @(Kn) (1+alpha(Kn).*Kn).*(1+4.*Kn./(1+Kn));
k_app_sample=k_app(Kn);
plot(Kn, k_app_sample(1,:),'b-', 'LineWidth', 2)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca,'FontSize',25);
xlim([1e-3 1e1]);
grid on;
xlabel('Knudsen number');
ylabel('Correcation factor Fapp');
title('Perm Correction factor (Florence et al, 2007)','FontSize',20);
%}

end