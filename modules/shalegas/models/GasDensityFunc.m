function [rho,rhoS,Z]=GasDensityFunc(method,parameters,varargin)
%{
Set up proper gas density function. rho(p)
Currently available:

1. Explicit method Z from explicit eqn (https://doi.org/10.1080/10916466.2016.1270305)
2. Implicit method Z from Peng-Robinson Equation
3. Table fitting. Fiting function with table of [P,rho]

One can check p-rho curve @ ...
http://www.unitrove.com/engineering/tools/gas/natural-gas-density

Arguments
---------
method -- 'Empirical' 'PR-EOS' 'Table'
Z      -- Gas compressibility factor, dimesnionless
M      -- average molecular weight, kg/mol
Temp   -- reservoir temperature, K
R      -- Boltzmann constant, J/(K.mol)
p_c    -- critical pressure of mix gas, Pa
T_c    -- critical temperature of mix gas, K
omega  -- PR-EOS parameters

Return
---------
[rho,rhoS] [func of rho(p), surface standard condition density of natural gas]
kg/m3

Author: Bin Wang(binwang.0213@gmail.com)
Date: Nov.2018
%}

opt = struct('TableData',[]);  % Compatibility only
opt = merge_options(opt, varargin{:});

MW=parameters.MW; %kg/mol
p_c=parameters.p_c; % pa
T_c=parameters.T_c;  %K
omega=parameters.PR_omega;
R=8.314;   %J/(K.mol)

Temp=parameters.T_Res;  %K 343.15
p_sc=parameters.p_sc;
T_sc=parameters.T_sc;
%T_sc=298.15;

if(strcmp(method,'Empirical'))
    %Reservoir Density
    T_pr  = @(T) T/T_c;
    p_pr  = @(p) p./p_c;
    Z     = @(p,T) 0.702*exp(-2.5*T_pr(T)).*p_pr(p).^2 ...
                 -5.524*exp(-2.5*T_pr(T)).*p_pr(p) ...
                 +0.044*T_pr(T)^2-0.164*T_pr(T)+1.15;
    
    %Surface Density
    Z_sc=Z(p_sc,T_sc);
end
    
if(strcmp(method,'PR-EOS'))
   %Analytical solution by book
   %Appendix B. Introductory Chemical Engineering Thermodynamics, Lira
   a = 0.457235*R^2 * T_c^2/p_c;
   b = 0.0777961*R*T_c/p_c; %Convert coefficient for MPa->Pa
   kappa = 0.37464 + 1.54226*omega - 0.26992*omega^2;
   alpha = @(T) (1 + kappa*(1 - (T/T_c).^0.5))^2;
   A = @(p,T) a*alpha(T).*p./(R^2)./(T^2);
   B = @(p,T) b.*p./R./T;
   
   %cubic function coefficient
   a2= @(p,T) B(p,T)-1; 
   a1= @(p,T) A(p,T)-3*B(p,T).^2-2.*B(p,T); 
   a0= @(p,T) -(A(p,T).*B(p,T)-B(p,T).^2-B(p,T).^3);
   cp= @(p,T) (3.*a1(p,T)-a2(p,T).^2)./3;
   cq= @(p,T) (2.*a2(p,T).^3-9.*a1(p,T).*a2(p,T)+27.*a0(p,T))./27;
   D=  @(p,T) 1/4.*cq(p,T).^2+1/27.*cp(p,T).^3;
   %%This is only valid for methane which D>0 only one root
   cP = @(p,T) (-cq(p,T)./2.0+D(p,T).^0.5).^(1/3.0);
   cQ = @(p,T) (-cq(p,T)./2.0-D(p,T).^0.5).^(1/3.0);
   Z  = @(p,T) (cP(p,T)+cQ(p,T))-a2(p,T)./3.0;
   
   %Surface Density
   Z_sc=Z(p_sc,T_sc);
end

if(strcmp(method,'Table'))%Z-factor from table by the 1st input parameters 
   P_Z_Table=opt.TableData;
   P_Z_Table = extendTab(P_Z_Table); % extend to constant values.
   Z=@(p,T) interpReg({P_Z_Table}, p, {':'});
   
   Z_sc=Z(p_sc,T_sc);
end

%{
p    = linspace(1*barsa, 400*barsa, 50);
figure('rend','painters','pos',[10 10 800 600]);
plot(convertTo(p, barsa), Z(p,Temp), 'bo-', 'LineWidth', 2)

T_pr  = @(T) T/T_c;
p_pr  = @(p) p./p_c;
Z1     = @(p,T) 0.702*exp(-2.5*T_pr(T)).*p_pr(p).^2 ...
             -5.524*exp(-2.5*T_pr(T)).*p_pr(p) ...
             +0.044*T_pr(T)^2-0.164*T_pr(T)+1.15;
hold on;
plot(convertTo(p, barsa), Z1(p,Temp), 'r^-', 'LineWidth', 2)

%CMG Data
fname=strcat(pwd,'\examples\data\benchmark_CMG_Fluid_Properties.csv');M = csvread(fname,1);
scatter(convertTo(M(:,1).*psia, barsa), M(:,2),80, 'ks');
hold off;
%plot(convertTo(p, psia), convertTo(mu(p),centi*poise),'bo-', 'LineWidth', 2)
set(gca,'FontSize',25);
xlabel('Pressure [Bar]')
ylabel('Z-factor')
grid on;
ylim([0.85 1]);
title('Gas Z-factor ','FontSize',20)
legend('Peng-Robinson','Empirical(Mahmoud, 2013)')
%}

%Gas density function
rho   = @(p) p.*MW./(Z(p,Temp).*R*Temp);
rhoS     = p_sc*MW/(Z_sc*R*T_sc);


end