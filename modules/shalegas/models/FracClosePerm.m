function [k_closefrac]=FracClosePerm(fracType,rigidity,p0,FracCellMask,parameters)
%{
Set up matrix perm reduction due to fracture clousure or proppant embeddment

Reference
---------
ARMA-2012-291 for propped fractures
Stiff log(Fcd)=-0.00011*S_close-0.0971
Meidum log(Fcd)=-0.00035*S_close+0.2396
Soft log(Fcd)=-0.00064*S_close-0.4585

SPE-184858-PA for unpropped fractures

Arguments
---------
Fcd   -- fracture conducivity kf*wf, md.ft
S_hmin -- horizontal minimum stress,Pa
S_hmax -- horizontal maximum stress,Pa
S_v    -- vertical overburden stress, Pa
p     -- fracture cell pressure,Pa
S_close -- fracture closure stress, Pa

Return
---------
[k_closefrac] permeability corrfection term due to fracture clousure

Author: Bin Wang(binwang.0213@gmail.com)
Date: Dec.2018
%}


S_v=38*mega*Pascal;   
S_hmin=parameters.S_hmin;  
S_hmax=parameters.S_hmax;   

if(strcmp(fracType,'HydraulicFrac'))
    S_close=@(p) convertTo(S_hmin, psia)-p;
else%Natural Frac
    S_close=@(p) convertTo((S_hmin+S_hmax)/2,psia)-p;
end

if(strcmp(fracType,'HydraulicFrac'))
    if(strcmp(rigidity,'Stiff'))
        Fcd0=10^(-0.00011.*S_close(convertTo(p0, psia))-0.0971);
        k_closefrac_func= @(p) 10.^(-0.00011.*S_close(p./psia)-0.0971)./Fcd0;
    elseif(strcmp(rigidity,'Medium'))
        Fcd0=10^(-0.00036.*S_close(convertTo(p0, psia))+0.2191);
        k_closefrac_func= @(p) 10.^(-0.00035.*S_close(p./psia)+0.2396)./Fcd0;
    elseif(strcmp(rigidity,'Soft'))
        Fcd0=10^(-0.0006.*S_close(convertTo(p0, psia))-0.4256);
        k_closefrac_func= @(p) 10.^(-0.00064.*S_close(p./psia)-0.4585)./Fcd0;
    end
else %Natural fractures
    if(strcmp(rigidity,'Stiff'))
        Fcd0=10.^(-0.793.*log(S_close(p0./psia))+4.5618);
        k_closefrac_func= @(p) 10.^(-0.793.*log(S_close(p./psia))+4.5618)./Fcd0;
    elseif(strcmp(rigidity,'Medium'))
        Fcd0=10.^(-0.89.*log(S_close(p0./psia))+5.0725);
        k_closefrac_func= @(p) 10.^(-0.89.*log(S_close(p./psia))+5.0725)./Fcd0;
    elseif(strcmp(rigidity,'Soft'))
        Fcd0=10.^(-1.041.*log(S_close(p0./psia))+6.0216);
        k_closefrac_func= @(p) 10.^(-1.041.*log(S_close(p./psia))+6.0216)./Fcd0;
    end
end

%Hydraulic fracture closure only occurs on the fracture 
k_closefrac= @(p) FracCellMask.*k_closefrac_func(p)+~FracCellMask;

%{
S_hmin=29*mega*Pascal;
S_hmax=34*mega*Pascal;
p    = linspace(34.5*barsa, p0, 20);
figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
%Stiff
S_close=@(p) convertTo(S_hmin, psia)-p;
Fcd0=10^(-0.0001*S_close(convertTo(p0, psia))-0.1082);
k_closefrac= @(p) 10.^(-0.0001*S_close(p./psia)-0.1082)./Fcd0;
plot(convertTo(p, barsa), k_closefrac(p),'bo-', 'LineWidth', 2)
hold on;
%Soft
S_close=@ (p) convertTo((S_hmin+S_hmax)/2,psia)-p;
Fcd0=10.^(-1.041.*log(S_close(p0./psia))+6.0216);
k_closefrac= @(p) 10.^(-1.041.*log(S_close(p./psia))+6.0216)./Fcd0;
plot(convertTo(p, barsa), k_closefrac(p),'r^-', 'LineWidth', 2)
hold off;
grid on;
set(gca,'FontSize',25);
xlabel('Pressure [Bar]')
ylabel('F_{app}')
xlim([34.5,p0/barsa]);
legend('Hydraulic Fracs','Natural Fracs')
title('Fracture Closure Correction Factor','FontSize',20)
%}

end