function fluid=setShaleGasFluid_Case2(G,rock,p0)
%{
General shale gas model manager
All of input parameters are distributed using parameters struct

Author: Bin Wang(binwang.0213@gmail.com)
Date: Dec.2018
%}

%% Input parameters table
parameters=struct();

%Unit converter
scf_ton=1/32*1e-6/1e-3; %[scf/ton->m3/kg]

%Basic parameters
parameters.C_Mat=4.4e-10; %[1/Pa] Matrix compressibility
parameters.T_Res=352;       %[K] Reservoir temperature
parameters.p_sc=1*atm;       %[Pa] Surface standrad pressure
parameters.T_sc=273.15 + 25; %[K] Surface standard temperature

%Fluid density and viscosity parameters
parameters.MW=16.04*1e-3;      %[kg/mol] Gas molecular Weight
parameters.p_c=46.4*barsa;     %[Pa] Gas Critical pressure
parameters.T_c=191;            %[K]  Gas Critical temperature
parameters.PR_omega=0.011;     %[-] PR-EOS parameters
parameters.rho_bulk=2500*kilogram/meter^3; %[kg/m3] matrix bulk density

%Adsorption parameters
parameters.p_langmuir=44.7*barsa;%[Pa] Langmuir pressure
parameters.v_langmuir=2.72e-3;   %[m3/kg] Langmuir volume

parameters.v_BET=124.53*scf_ton; %[m3/kg] BET volume
parameters.C_BET=36.63;%[-] BET C constant 
parameters.n_BET=4.03;   %[-] BET n constant

%Geomechanics SPE-170830-PA
parameters.m=0.5;                %[-] Gangi coefficient
parameters.Pc=103*mega*Pascal;   %[MPa] overburden vertical stress
parameters.P1=180*mega*Pascal;   %[MPa] fracture closure stress
parameters.alpha=0.5;            %[-] Biot constant

%Fracture closure and embeddment
parameters.S_hmin=29*mega*Pascal;   %[MPa] min horizontal stress
parameters.S_hmax=34*mega*Pascal;   %[MPa] max horizontal stress


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
[G,deltaM,deltaHF,deltaNF]=setDomainDeltaFunc(G);

%% [Optional] Shale gas adsorption term in the matrix
fluid.rho_bulk= parameters.rho_bulk;%kg/m3
parameters.rhoGS=fluid.rhoGS;
fluid.mG_ad=MatrixAdsorptionFunc('Langmuir',G.FracCellMask,parameters);

%% [Optional] Shale gas appraent perm for gas slippage flow in the matrix
%fluid.kG_app=MatrixApparentPerm('Civan',rock,fluid.muG,G.FracCellMask,parameters);

%% [Optional] Micro-fracture closure for geomechanics effect in the matrix
%fluid.k_gangi=MatrixGangiPerm(G.FracCellMask);

%% [Optional] Hydraulic fracture and natural fracture closure for geomechanics effect
fluid.k_hydraulicfrac=FracClosePerm('HydraulicFrac','Stiff',p0,G.HFCellMask,parameters);
fluid.k_naturalfrac=FracClosePerm('NaturalFrac','Soft',p0,G.NFCellMask,parameters);

%% [Optional] Non-darcy Forchheimer in the fracture
%fluid.k_nondarcy=FracNonDarcy(rock,fluid.rhoG,fluid.muG,G.FracCellMask);

end