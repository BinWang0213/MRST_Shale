function viscMix = viscLBC(c, MW, TinK, Pci, Tci, Vci, zi)
%VISCLBC Summary of this function goes here
%   Detailed explanation goes here
    %convert Pci from Pa to atm
    Pci = Pci ./ 101325.0;      %atm

    Tri  = TinK./Tci;
    eta_i = (Tci./((MW.^3.0).*(Pci.^4.0))).^(1.0/6.0);
    
    %testing if Tri<1.5 is equivalent to checking if mu_i*etai < 4.9774e-04
    mu0_iEta_i = 34.0e-5 .* (Tri.^0.94);
    mu0_iEta_i(mu0_iEta_i>4.97742472e-04) = ...
        17.78e-5 .* ((4.58 .*Tri(Tri>1.5) - 1.67).^(5.0/8.0));
    
    mu0_i = mu0_iEta_i ./ eta_i;
    
    %Herning and Zipperer equation
    sqrtMW = sqrt(MW);
    mu0_mix = ((zi .* mu0_i)*sqrtMW') / (zi*sqrtMW');
    
    rho_c = 1.0 ./ (zi*Vci'); % pseudo-critical density
    rho_r = c / rho_c;        % pseudo-reduced density. Both num & denom are in mol/m^3  
    
    etaMix = ((zi*Tci')/(((zi*MW')^3.0)*((zi*Pci')^4.0)))^(1.0/6.0);
    viscMix = mu0_mix + ((((((0.0093324*rho_r-0.040758)*rho_r+0.058533) ...
              *rho_r+0.023364)*rho_r+0.1023)^4.0)-1.0e-4)/etaMix;
    viscMix = viscMix * 0.001; %convert viscosity from cp to Pa.s
end

