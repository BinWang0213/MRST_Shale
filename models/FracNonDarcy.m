function [k_nondarcy]=FracNonDarcy(rock,rho,mu,FracCellMask)
%{
Set up non-darcy flow for hydarulic fractures

Reference
---------
SPE-88536-MS
SPE-163609-PA
SPE-132093-MS

k_nondarcy=beta*k/miu^2
B_face=k_nondarcy*Trans*grad(p)
F_ND_face=-1+sqrt(1-4B)/(2B)

Arguments
---------
beta   -- Forchheimer parameter, 1/m
F_ND   -- Non-Darcy flow resistance factor, <1
B      -- parameter to calculate F_ND, >0
Trans  -- Cell face transmisbility 

Return
---------
[k_closefrac] permeability corrfection term due to fracture clousure

Author: Bin Wang(binwang.0213@gmail.com)
Date: Dec.2018
%}

%Forchheimer parameter, perm unit = md
beta=3.2808*1.485e9./(rock.perm/(milli*darcy)).^1.021;

%Fnd non-darcy correction coefficient
k_nondarcy_func= @(p) rho(p).*beta.*rock.perm./mu(p).^2;

%k_nondarcy= @(p) FracCellMask.*k_nondarcy_func(p)+~FracCellMask;
k_nondarcy= @(p) k_nondarcy_func(p);


%{
%Verified with Evans RD, Civan F. Characterization of non-Darcy multiphase flow
% Page B-9 Fig 4.5 beta=1/ft K=md
perm = linspace(1e-2*milli*darcy, 1e6*milli*darcy, 20);
figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
beta=3.2808*1.485e9./(perm/(milli*darcy)).^1.021;
plot(convertTo(perm,milli*darcy), beta.*1/3.2808,'bo-', 'LineWidth', 2)
set(gca,'FontSize',25);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('K[md]')
ylabel('\beta (ft^{-1})')
title('Correlation of \beta vs K','FontSize',20)

%F_NB
B=linspace(1e-15,1,20);
F_ND1=@(B) (-1+(1+4.*B).^0.5)./(2.*B);
F_ND2=@(B) 2./(1+(1+4.*B).^0.5);

figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
set(gca,'FontSize',25);
plot(B, F_ND1(B),'bo-', 'LineWidth', 2)
hold on;
plot(B, F_ND2(B),'r^-', 'LineWidth', 2)
hold off;
xlabel('B')
ylabel('F_ND')
title('Non-Darcy flow resistance factor','FontSize',20)
%}

end