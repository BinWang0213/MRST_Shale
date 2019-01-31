function fluid=setShaleGasFluid_Case1(G,rock,p0)
%{
General shale gas model manager
All of input parameters are distributed into each sub-function

Author: Bin Wang(binwang.0213@gmail.com)
Date: Dec.2018
%}

%% Input parameters table
parameters=struct();

%Unit converter
lb_ft3=16.0184634;
scf_ton=1/32*1e-6/1e-3; %[scf/ton->m3/kg]
gmoles_lb=5.97696667e-4; %[scf/ton->gmoles/lb]

%Basic parameters
parameters.C_Mat=1e-6*psia^-1; %[1/Pa] Matrix compressibility
parameters.T_Res=(130 - 32)*5/9 + 273.15; %[K] Reservoir temperature
parameters.p_sc=1*atm;       %[Pa] Surface standrad pressure
parameters.T_sc=273.15 + 15; %[K] Surface standard temperature

%Fluid density and viscosity parameters
parameters.MW=16.04*1e-3;      %[kg/mol] Gas molecular Weight
parameters.p_c=45.4*atm;     %[Pa] Gas Critical pressure
parameters.T_c=190.6;            %[K]  Gas Critical temperature
parameters.rho_c=162.66;        %[kg/m3] Gas Critical density
parameters.PR_omega=0.01142;     %[-] PR-EOS parameters

%Adsorption parameters
parameters.p_langmuir=1300*psia;%[Pa] Langmuir pressure
parameters.v_langmuir=130*scf_ton;  %[m3/kg] Langmuir volume
parameters.rho_bulk=156.075*lb_ft3; %[kg/m3] matrix bulk density


%% [Necessary] Black-oil fluid properties
fluid = initSimpleADIFluid('phases','G', ...
                           'cR',parameters.C_Mat, ...
                           'pRef',0*psia); 

%% [Necessary] Fluid model
%fname=strcat(pwd,'\examples\data\benchmark_CMG_Fluid_Properties.csv');
fname='CMG_Fluid_Properties.csv';
M = csvread(fname,1); % P,Z-factor,viscosity
P_Z=[M(:,1).*psia M(:,2).*1];
P_mu=[M(:,1).*psia M(:,3).*centi*poise];

[fluid.rhoG,fluid.rhoGS]=GasDensityFunc('Table',parameters,'TableData',P_Z);
%[fluid.rhoG,fluid.rhoGS]=GasDensityFunc('Empirical',parameters);
%[fluid.rhoG,fluid.rhoGS]=GasDensityFunc('PR-EOS',parameters);
fluid.bG=@(p) fluid.rhoG(p)./fluid.rhoGS;
fluid.rhoGS=fluid.rhoGS;
%fluid.muG=GasViscFunc('Lee',fluid.rhoG,parameters);
%fluid.muG=GasViscFunc('Jossi',fluid.rhoG,parameters);
fluid.muG=GasViscFunc('Table',fluid.rhoG,parameters,'TableData',P_mu);


%% Genetrate Frac Cell mask to switch physics for matrix and fracture
[deltaM,deltaHF,deltaNF]=setDomainDeltaFunc(G);

%% [Optional] Shale gas adsorption term in the matrix
G.FracCellMask(:)=0; %CMG has adsorptin for all grids

fluid.rho_bulk = parameters.rho_bulk;%kg/m3
parameters.rhoGS=fluid.rhoGS;
fluid.mG_ad=MatrixAdsorptionFunc('Langmuir',G.FracCellMask,parameters);

%% [Optional] Shale gas appraent perm for gas slippage flow in the matrix
%fluid.kG_app=MatrixApparentPerm('Civan',rock,fluid.muG,G.FracCellMask,parameters);

%% [Optional] Micro-fracture closure for geomechanics effect in the matrix
%fluid.k_gangi=MatrixGangiPerm(G.FracCellMask);

%% [Optional] Hydraulic fracture and natural fracture closure for geomechanics effect
%fluid.k_hydraulicfrac=FracClosePerm('HydraulicFrac','Stiff',p0,G.FracCellMask);
%fluid.k_naturalfrac=FracClosePerm('NaturalFrac','Soft',p0,G.FracCellMask);

%% [Optional] Non-darcy Forchheimer in the fracture
%fluid.k_nondarcy=FracNonDarcy(rock,fluid.rhoG,fluid.muG,G.FracCellMask);

end