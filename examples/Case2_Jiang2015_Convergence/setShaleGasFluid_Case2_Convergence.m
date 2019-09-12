function fluid=setShaleGasFluid_Case2_Convergence(G,rock,p0)
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
parameters.C_Mat=0.5e-4/barsa; %[1/Pa] Matrix compressibility
parameters.T_Res=343.15;       %[K] Reservoir temperature
parameters.p_sc=1*atm;       %[Pa] Surface standrad pressure
parameters.T_sc=273.15 + 25; %[K] Surface standard temperature

%Fluid density and viscosity parameters
parameters.MW=16.04*1e-3;      %[kg/mol] Gas molecular Weight
parameters.p_c=46.4*barsa;     %[Pa] Gas Critical pressure
parameters.T_c=191;            %[K]  Gas Critical temperature
parameters.PR_omega=0.011;     %[-] PR-EOS parameters
parameters.rho_bulk=2500*kilogram/meter^3; %[kg/m3] matrix bulk density

%Adsorption parameters
parameters.p_langmuir=40*barsa;%[Pa] Langmuir pressure
parameters.v_langmuir=0.018;   %[m3/kg] Langmuir volume


%% [Necessary] Black-oil fluid properties
fluid = initSimpleADIFluid('phases','G', ...
                           'cR',parameters.C_Mat, ...
                           'pRef',0*barsa); 

%% [Necessary] Fluid model
%[fluid.rhoG,fluid.rhoGS]=GasDensityFunc('Empirical',parameters);
[fluid.rhoG,fluid.rhoGS]=GasDensityFunc('PR-EOS',parameters);
fluid.bG=@(p) fluid.rhoG(p)./fluid.rhoGS;
fluid.rhoGS=fluid.rhoGS;
fluid.muG=GasViscFunc('Lee',fluid.rhoG,parameters);

%% Genetrate Frac Cell mask to switch physics for matrix and fracture
%[G,deltaM,deltaHF,deltaNF]=setDomainDeltaFunc(G);

%% [Optional] Shale gas adsorption term in the matrix
fluid.rho_bulk= parameters.rho_bulk;%kg/m3
parameters.rhoGS=fluid.rhoGS;
%fluid.mG_ad=MatrixAdsorptionFunc('Langmuir',G.FracCellMask,parameters);

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