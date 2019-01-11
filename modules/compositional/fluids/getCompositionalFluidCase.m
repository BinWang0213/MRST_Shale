function [fluid, info] = getCompositionalFluidCase(name, varargin)
% Grab bag of different fluids for use in examples etc
    switch(lower(name))
        case 'lumped_1'
            % From SPE 79691 ex 5
            T_c = [189.515, 304.200, 387.607, 597.497, 698.515, 875.00];
            P_c = [45.2012, 72.90, 40.4196, 33.0150, 17.4525, 11.5372]*atm;
            mw = [16.1594, 44.01, 45.5725, 117.740, 248.827, 481.520]/1000;
            acc = [0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141];
            Z_c = [0.28981, 0.27055, 0.27588, 0.25668, 0.21967, 0.18250];
            
            V_c = Z_c.*8.314.*T_c./P_c;
            
            names = {'N2/CH4', 'CO2', 'C2-5', 'C6-13', 'C14-24', 'C25-80'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            bic = [0.11883, 0, 0, 0,0, 0;...
                   0.00070981, 0.15,0, 0, 0, 0; ...
                   0.00077754, 0.15,0, 0, 0, 0; ...
                   0.01, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0];
               
            bic = bic + tril(bic, -1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 1, 0, 0, 0, 0], ...
                            'initial',   [0.463, 0.01640, 0.20520, 0.19108, 0.08113, 0.04319], ...
                            'pressure', 225*atm, ...
                            'temp', 387.45);
        case 'spe5'
            % Fifth SPE benchmark
            T_c = [343, 665.7, 913.4, 1111.8, 1270.0, 1380.0]*Rankine;
            P_c = [667.8, 616.3, 436.9, 304.0, 200.0, 162.0]*psia;
            mw = [16.040, 44.100, 86.180, 142.290, 206.000, 282.000]/1000;
            acc = [0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500];
            Z_c = [0.290, 0.277, 0.264, 0.257, 0.245, 0.235];
            
            V_c = Z_c.*8.314.*T_c./P_c;
            
            names = {'C1', 'C3', 'C6', 'C10', 'C15', 'C20'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(5, 1) = 0.05;
            bic(6, 1) = 0.05;
            bic(5, 2) = 0.005;
            bic(6, 2) = 0.005;
            
            bic = bic + tril(bic, -1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0.77, 0.20, 0.03, 0, 0, 0], ...
                            'initial',   [0.5, 0.03, 0.07, 0.20, 0.15, 0.05], ...
                            'pressure', 4000*psia, ...
                            'temp', 344.26);

        case 'mix1'
            % 5 component light-gas mixture
            names = {'Methane', 'Ethane', 'n-Propane', 'Nitrogen', 'CarbonDioxide'};
            info = makeInfo('injection', [0, 0, 0, 0, 1.0], ...
                            'initial',   [0.982, 0.0066, 0.0001, 0.0091, 0.0022], ...
                            'pressure', 4000*psia, ...
                            'temp', 353.15); 
            fluid = TableCompositionalFluid(names);
            
        case 'mix1_3comps'
            %OMO: Edit added a couple of hydrocarbon mixtures
            % 5 component light-gas mixture
            names = {'Methane', 'Ethane', 'n-Propane'};
%             isotherm = [1057.9, 322.1, 0; 24.18, 4.0187, -0.0576128; 21.235, -0.79523, 0.020037]'*mega*Pascal;
            isotherm = [1.6926E+01,5.1508E-06,0; 7.3126E-01,1.2042E-07,-1.7197E-15; 9.3707E-01,-3.5165E-08,8.8618E-16]';
            fluid = NanoCompositionalFluid(names,isotherm);    
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 2) = 0.005;
            bic(1, 3) = 0.01;
            bic(2, 3) = 0.005;
            
            bic = bic + triu(bic, 1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 0, 1.0], ...
                            'initial',   [0.991,0.0088,0.0002], ...
                            'pressure', 5466.482220000024*psia, ...
                            'temp', 355.84); 
            
        case 'mix2'
            % 5 component intermediate-gas mixture
            names = {'Methane', 'Ethane', 'n-Propane', 'n-Butane', 'n-Pentane'};
            info = makeInfo('injection', [0, 0, 0, 0, 1.0], ...
                            'initial',   [0.749, 0.097, 0.086, 0.048,0.02], ...
                            'pressure', 4000*psia, ...
                            'temp', 353.15); 
            fluid = TableCompositionalFluid(names);
            
        case 'mix3'
            % 5 component heavy gas mixture
            names = {'Methane', 'Ethane', 'n-Propane', 'n-Butane', 'n-Pentane'};
            info = makeInfo('injection', [0, 0, 0, 0, 1.0], ...
                            'initial',   [0.538, 0.164, 0.127, 0.105,0.066], ...
                            'pressure', 4000*psia, ...
                            'temp', 353.15); 
            fluid = TableCompositionalFluid(names);            
  
        case 'mix3co2'
            % Mix3 + CO2
            names = {'Methane', 'Ethane', 'n-Propane', 'n-Butane', 'n-Pentane', 'CarbonDioxide'};
            info = makeInfo('injection', [0.269, 0.081, 0.063, 0.053, 0.034, 0.500], ...
                            'initial',   [0.54, 0.166, 0.127, 0.107,0.059, 0.001], ...
                            'pressure', 4000*psia, ...
                            'temp', 353.15); 
            fluid = TableCompositionalFluid(names);
            
        case 'condensate'
            % Condensate at 4000 psia and 165 F
            names = {'Methane', 'Ethane', 'n-Butane', 'n-Heptane', 'n-Tetradecane'};
            info = makeInfo('injection', [0, 0, 0, 0, 1.0], ...
                            'initial',   [0.6442, 0.1723, 0.1322, 0.0372, 0.0141], ...
                            'pressure', 4000*psia, ...
                            'temp', 347.039); 
            fluid = TableCompositionalFluid(names);    
            
        case 'drygas'
            % dry at 8000 psia and 210 F
            names = {'Methane', 'Ethane', 'n-Butane', 'n-Heptane', 'n-Tetradecane'};
            info = makeInfo('injection', [0, 0, 0, 0, 1.0], ...
                            'initial',   [0.6442, 0.1723, 0.1322, 0.0372, 0.0141], ...
                            'pressure', 8000*psia, ...
                            'temp', 372.039); 
            fluid = TableCompositionalFluid(names);              

        case 'marat'
            % dry at 4000 psia and 210 F
            names = {'Methane', 'Ethane', 'n-Butane', 'n-Heptane', 'n-Decane'};
            info = makeInfo('injection', [0, 0, 0, 0, 1.0], ...
                            'initial',   [0.6442, 0.1723, 0.1322, 0.0372, 0.0141], ...
                            'pressure', 4000*psia, ...
                            'temp', 372.039); 
            fluid = TableCompositionalFluid(names); 
                        
                        
        case 'fourmix'
            names = {'Methane', 'Nitrogen', 'n-Pentane', 'n-Decane'};
            
            info = makeInfo('injection', [0, 1, 0, 0], ...
                            'initial',   [0.1, 0.0, 0.45, 0.45], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 50);
                        
            fluid = TableCompositionalFluid(names);
            
        case 'bakken'
            % Bakken shale composition from Nojabaei dissertation 2015
            T_c = [335.336	549.969	665.97	759.208	875.479	1053.25	1332.095	1844.491].*Rankine();
            P_c = [655.02	721.99	615.76	546.46	461.29	363.34	249.61	190.12]*psia;
            mw = [16.535	30.433	44.097	58.124	78.295	120.562	220.716	443.518]/1000;
            acc = [0.0102	0.1028	0.152	0.1894	0.2684	0.4291	0.7203	1.0159];
%             Z_c = [0.290, 0.277, 0.264, 0.257, 0.245, 0.235];

            
            V_c = [1.58	2.34	3.25	4.11	5.39	8.81	15.19	36].*0.06242795996802./1000;
            Z_c = P_c.* V_c./8.314./T_c;
            Parachor = [74.8	107.7	151.9	189.6	250.2	350.2	590	1216.8];
            
            names = {'C1','C2','C3','C4','C5-C6','C7-C12','C13-C21','C22-C80'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(1, 1:8) = [0	0.005	0.0035	0.0035	0.0037	0.0033	0.0033	0.0033];
            bic(2, 1:8) = [0	0	0.0031	0.0031	0.0031	0.0026	0.0026	0.0026];

            bic = bic + triu(bic, 1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [1	0 0 0 0 0 0 0 ], ...
                            'initial',   [0.36736	0.14885	0.09334	0.05751	0.06406	0.15854	0.0733	0.03704], ...
                            'pressure', 6840*psia, ...
                            'temp', 388.706);
            
        case '20components'
            names = {'Methane', 'Nitrogen', 'n-Pentane', 'n-Decane', 'CycloHexane', ...
                     'CarbonMonoxide', '1-Butene', 'Ammonia', 'HydrogenSulfide', 'IsoButane', ...
                     'n-Dodecane', 'n-Heptane', 'n-Hexane', 'n-Nonane', 'n-Octane', ...
                     'p-Xylene', 'n-Propane', 'n-Undecane', 'Oxygen', 'Water'};
            info = makeInfo('injection', (1:20)/20, ...
                            'initial',   ones(1, 20)/20, ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 50);
                        
            fluid = TableCompositionalFluid(names);

        case {'liquid_initial', 'vapor_initial'}
            names = {'Oxygen', 'Water'};
            if strcmpi(name, 'liquid_initial')
                inj =  [1, 0];
                init = [0, 1];
            else
                inj =  [0, 1];
                init = [1, 0];
            end            
            info = makeInfo('injection', inj, ...
                            'initial',   init, ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 30);
                        
            fluid = TableCompositionalFluid(names);
            fluid.names{2} = 'H2O';
            fluid.names{1} = 'O2';
        case 'watertracer'
            fluid = TableCompositionalFluid({'water', 'water'});
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0, 1], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 30);
            fluid.names = {'Tracer1', 'Tracer2'};
        case 'simple'
            names = {'Methane', 'CarbonDioxide', 'n-Decane'};
            
            info = makeInfo('injection', [0.1, 0.9, 0], ...
                            'initial',   [0.3, 0.1, 0.6], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalFluid(names(order));
        case 'verysimple'
            names = {'CarbonDioxide', 'n-Decane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.1, 0.9], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
        case 'onlydecane'
            names = {'n-Decane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1], ...
                            'initial',   [1], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
        case 'onlymethane'
            names = {'Methane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1], ...
                            'initial',   [1], ...
                            'pressure', 5000*psia, ...
                            'temp', 366.48);
        case 'c1c2'
            names = {'Methane','Ethane'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [0.999 0.001], ...
                            'initial',   [0.999 0.001], ...
                            'pressure', 5000*psia, ...
                            'temp', 366.48); 
                        
        case 'watermethane'
            names = {'Methane', 'Water'};
            fluid = TableCompositionalFluid(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.01, 0.99], ...
                            'pressure', 5000*psia, ...
                            'temp', 366.48);
        otherwise
            error('Unknown case')
    end
end

function info = makeInfo(varargin)
    info = struct('injection', [], 'initial', [],...
                  'pressure',  [], 'temp',    []);
    info = merge_options(info, varargin{:});
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
