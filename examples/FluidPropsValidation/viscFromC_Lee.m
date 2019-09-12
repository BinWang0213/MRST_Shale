function visc = viscFromC_Lee(c)
%Computes viscosity from density using Lee et al., 1966
%   Paper is titled- The Viscosity of Natural Gases 
%               
    global molWeight;
    global TinK;
    TinR = (TinK-273.15)*1.8 + 491.67;
    MW_g = molWeight*1000.0;
    dens = c*molWeight/1000.0; %converting to g/cc
    X = 2.57 + 1914.5/TinR + 0.0095*MW_g;
    Y = 1.11 + 0.04*X;
    K = 1e-7*(7.77+0.0063*MW_g)*(TinR^1.5)/(122.4+12.9*MW_g+TinR);
    visc = K*exp(X*(dens^Y)); %in Pa.s
end