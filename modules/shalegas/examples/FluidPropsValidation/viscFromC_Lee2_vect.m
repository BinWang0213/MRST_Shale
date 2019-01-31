function visc = viscFromC_Lee2_vect(c, MWi, TinK, varargin)
%Computes viscosity from density using Lee et al., 1966
%   Paper is titled- The Viscosity of Natural Gases 
%      c, MWi, TinK, Pci, Tci, Vci, zi         
    
    if (nargin == 3)
        MW_kg = MWi/1000.0; %molecular weight in kg/mol
    else
        xi = varargin{1};
        MW_kg = (xi*MWi')./1000.0; %apparent molecular weight in kg/mol    
    end

    TinR = (TinK-273.15)*1.8 + 491.67;
    MW_g = MW_kg*1000.0;   %apparent molecular weight in g/mol
    
    dens = c.*MW_kg/1000.0; %converting to g/cc
    X = 3.448 + 986.4/TinR + 0.01009*MW_g;
    Y = 2.447 - 0.2224*X;
    K = 1e-7*(9.379+0.01607*MW_g)*(TinR^1.5)/(209.2+19.26*MW_g+TinR);
    visc = K.*exp(X*(dens.^Y)); %in Pa.s
end