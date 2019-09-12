function [m_ad]=MatrixAdsorptionFunc(method,FracCellMask,parameters)
%{
Set up adsorption accumation func
Currently available:
1. Langmuir ()
2. BET (SPE-170801-PA)

Arguments
---------
method -- 'Langmuir' 'BET'
PL     -- Langmuir pressure of gas mixture, pa
VL     -- Langmuir volume of gas mixture, m3/kg

Ps     -- pseudo saturation pressure, MPa


Return
---------
[m_ad] [func of accumation term considering adsorption]

Author: Bin Wang(binwang.0213@gmail.com)
Date: Nov.2018
%}

rhoS=parameters.rho_bulk;
rhoG_SC=parameters.rhoGS;

if(strcmp(method,'Langmuir'))
    PL=parameters.p_langmuir;% Pa
    VL=parameters.v_langmuir; %m3/kg

    m_ad_func = @(p) p.*VL./(PL+p);
    %{
    p    = linspace(1*barsa, 300*barsa, 50);
    figure('rend','painters','pos',[10 10 800 600]);
    plot(convertTo(p, barsa), m_ad_func(p),'bo-', 'LineWidth', 2)
    set(gca,'FontSize',25);
    xlabel('Pressure [Bar]')
    ylabel('Adsorbed Gas [m^3/kg]')
    title('Langmuir Isotherm','FontSize',20)
    %}
end

if(strcmp(method,'BET')) 
    T=parameters.T_Res;
    C=parameters.C_BET;
    vm=parameters.v_BET;
    n=parameters.n_BET;
    
    
    p_s= exp(7.7437-1306.5485/(19.4362+T))*mega*Pascal;
    p_r= @(p) p./p_s;
    m_ad_func= @(p) vm*C.*p_r(p)./(1-p_r(p))...
                   .*(1-(n+1).*p_r(p).^n+n.*p_r(p).^(n+1))...
                   ./(1+(C-1).*p_r(p)-C.*p_r(p).^(n+1));
    %{ 
    %Match plot with sample 1 in Fig.3c of SPE-170801-PA
    T= 327.594;
    scf_ton=1/32*1e-6/1e-3;
    vm=49.01*scf_ton;
    C=24.56;
    n=4.46;
    p    = linspace(1*barsa, 300*barsa, 50);
    p_s= exp(7.7437-1306.5485/(19.4362+T))*mega*Pascal;
    p_r= @(p) p./p_s;
    m_ad_func= @(p) vm*C.*p_r(p)./(1-p_r(p))...
                   .*(1-(n+1).*p_r(p).^n+n.*p_r(p).^(n+1))...
                   ./(1+(C-1).*p_r(p)-C.*p_r(p).^(n+1));
               
    VL=100.6*scf_ton;
    PL=1144*psia;
    m_ad_func_langmuir = @(p) p.*VL./(PL+p);

    figure('rend','painters','pos',[10 10 800 600]);
    plot(convertTo(p, barsa), m_ad_func(p)/1,'bo-', 'LineWidth', 2)
    hold on;
    plot(convertTo(p, barsa), m_ad_func_langmuir(p)/1,'r^-', 'LineWidth', 2)
    hold off;
    set(gca,'FontSize',25);
    xlabel('Pressure [Bar]')
    ylabel('Adsorbed Gas [m^3/kg]')
    grid on;
    title('Langmuir/BET Isotherm','FontSize',20)
    legend('BET','Langmuir')
    %}
end

%Adsorption only occurs on the matrix
m_ad = @(p) ~FracCellMask.*rhoS.*rhoG_SC.*m_ad_func(p);

end