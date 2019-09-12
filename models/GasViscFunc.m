function [mu]=GasViscFunc(method,rho,parameters,varargin)
%{
Set up proper gas viscosity function. miu(p)


1. Lee 1966, SPE-1340-PA
sample calculation
https://koyapete.weebly.com/uploads/1/3/6/4/13645543/fundamentals_of_petroleum_lec._11.pdf
verificatoin curve
https://checalc.com/solved/gasVisc.html

Application range:
6.89 bar < p < 551.5806 bar, 
310.928 < T(K) < 444.261 
0.90 < CO2 (mol%) < 3.20 
0.0 < N2 (mol%) < 4.80

2. Jossi
https://www.e-education.psu.edu/png520/m19_p4.html
P17 1.2.2.1 @ Viscosity Study of Hydrcarbon Fluids at Reservoir Conditions Modeling and
Measurements

Arguments
---------
rho    -- func of natural density whic depends on pressure, g/cm3
M      -- average molecular weight, kg/mol
T      -- reservoir temperature, Rankine

Return
---------
[fe]:body force vector (size=nx1)

Author: Bin Wang(binwang.0213@gmail.com)
Date: Nov.2018
%}

opt = struct('TableData',[]);  % Compatibility only
opt = merge_options(opt, varargin{:});

if(strcmp(method,'Lee'))
    M=parameters.MW*1e3; %kg/mol->g/mol
    T=parameters.T_Res*1.8;  %K->Rankine
    
    K=(9.379+0.01607*M)*(T^1.5)/(209.2+19.26*M+T);
    X=3.448+986.4/T+0.01009*M;
    Y=2.447-0.2224*X;
    
    mu = @(p) 1e-4*K.*exp(X.*(1e-3.*rho(p)).^Y)*centi*poise; %[cp]
end
if(strcmp(method,'Jossi'))
    p_c=convertTo(parameters.p_c,atm); % pa
    T_c=parameters.T_c;  %K
    MW=parameters.MW*1e3; %kg/mol->g/mol
    rho_c=parameters.rho_c; % kg/m3
    T=parameters.T_Res; %K
    
    xi=T_c^(1/6)/(MW^(1/2)*p_c^(2/3));
    rho_r= @(p) rho(p)./rho_c;
    
    %Using Lee's model to calculate mu0
    mu0=1e-4*(17.94+0.0321*MW)*T^(1.5)/(1.8*T+75.4+13.9*MW)*0.84;
    
    mu=@(p) ( mu0 - 1e-4...
        +1/xi.*(0.1023+0.023364.*rho_r(p)...
        +0.058533.*rho_r(p).^2 ...
        -0.040758.*rho_r(p).^3 ...
        +0.0093724.*rho_r(p).^4).^4 )*centi*poise; %[cp]
end

if(strcmp(method,'Table')) %Fist parameter is a [Nx2] Table
    P_mu_Table=opt.TableData;
    P_mu_Table = extendTab(P_mu_Table); % extend to constant values.
    mu=@(p) interpReg({P_mu_Table}, p, {':'});
end

%{
    p    = linspace(1*barsa, 400*barsa, 50);
    figure('rend','painters','pos',[10 10 800 600]);
    plot(convertTo(p, barsa), convertTo(mu(p),centi*poise),'bo-', 'LineWidth', 2)
    %plot(convertTo(p, psia), convertTo(mu(p),centi*poise),'bo-', 'LineWidth', 2)
    hold on;
    %CMG Data
    fname=strcat(pwd,'\examples\data\benchmark_CMG_Fluid_Properties.csv');M = csvread(fname,1);
    scatter(convertTo(M(:,1).*psia, barsa), M(:,3),80, 'ks');
    hold off;
    set(gca,'FontSize',25);
    xlabel('Pressure [Bar]')
    ylabel('Viscosity [cp]')
    grid on;
    title('Gas visocity (Lee, 1966)','FontSize',20)
%}


end