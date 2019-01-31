clc
clear
close all

mrstModule add compositional deckformat ad-core ad-props 

G =  cartGrid([1 1],[1 1]);
G = computeGeometry(G);
rock = makeRock(G, 50*milli*darcy, 0.25);

casename = 'onlymethane';
[fluid, info] = getCompositionalFluidCase(casename);
% p = info.pressure
T = info.temp
z = 1
load('CMGFluidProperties.mat')

p = CMGFluidProperties.Ppsi .* psia;
Z_g_CMG = CMGFluidProperties.Zfactor; 
visc_CMG = CMGFluidProperties.viscositycp; 
% T = zeros(numel(p),1).*info.temp;
% z = zeros(numel(p),1).*1;
eosname = 'pr';  %'srk','rk','prcorr'
EOSModel = EquationOfStateModel(G, fluid, eosname);
[L, x, y, Z_L, Z_V,mu_MRST] = deal(zeros(numel(p),1));

isLiquid = L;
nkr = 2;
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);
model = OverallCompositionCompositionalModel(G, rock, flowfluid, fluid, 'water', false);

for i=1:numel(p)
    [L(i), x(i), y(i), Z_L(i), Z_V(i)] = standaloneFlash(p(i), T, z, EOSModel);
    mu_MRST(i) = 1000.*model.PropertyModel.computeViscosity(p(i), y(i), Z_V(i), T, isLiquid(i) );
end


eosname = 'prcorr';  %'srk','rk','prcorr'
EOSModel = EquationOfStateModel(G, fluid, eosname);
[L2, x2, y2, Z_L2, Z_V2] = deal(zeros(numel(p),1));

for i=1:numel(p)
    [L2(i), x2(i), y2(i), Z_L2(i), Z_V2(i)] = standaloneFlash(p(i), T, z, EOSModel);
end

Pci = 4604000;
Tci = 190.5800;  
accFact_i = 0.0108; 
MWi = 16.043;
Vci = 99.3e-6; %m^3/mol
RR = 8.3145;

Zg_OMO = zFactPureAllCVs(p,T,Pci,Tci,accFact_i,RR,numel(p));




c= p./(Zg_OMO*RR*T);

[viscLee, viscLBC1] = deal(zeros(numel(p),1));
for i=1:numel(p)
    viscLee(i)  = 1000.*viscFromC_Lee2(c(i), MWi, T);
    viscLBC1(i) = 1000.*viscLBC(c(i), MWi, T, Pci, Tci, Vci, 1);
end




p = p./psia;

figure(1)
plot(p, Z_g_CMG,'r*',p, Zg_OMO,'g', p, Z_V, 'b--')
% plot(p, Z_V2,'r', p, Z_V, 'b--')
title('Comparison of Z-factor Models')
xlabel('Pressure, psi')
ylabel('Z-factor')
legend({'Z_g CMG','Z_g OMO','Z_g MRST'})

figure(2)
plot(p, visc_CMG,'r*',p, viscLBC1,'g', p, viscLee, 'b--', p, mu_MRST)
title('Comparison of Viscosity Models')
xlabel('Pressure, psi')
ylabel('Viscosity, cp')
legend({'visc CMG','visc LBC','visc Lee','visc MRST'})