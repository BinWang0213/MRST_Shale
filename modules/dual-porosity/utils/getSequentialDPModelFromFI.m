function model = getSequentialDPModelFromFI(fimodel, varargin)
% For a given fully implicit model, output the corresponding pressure/transport model 
    if isa(fimodel, 'SequentialPressureTransportDPModel')
        % User gave us a sequential model! We do not know why, but in that
        % case we just return it straight back.
        model = fimodel;
        return
    end
    rock  = fimodel.rock;
    fluid = fimodel.fluid;
    rock_matrix  = fimodel.rock_matrix;
    fluid_matrix = fimodel.fluid_matrix;
    dp_info = fimodel.dp_info;
    
    G     = fimodel.G;
    
    switch lower(class(fimodel))
        case 'twophaseoilwaterdpmodel'
            pressureModel  = PressureOilWaterDPModel(G, rock, fluid, ...
                                                    rock_matrix, fluid_matrix, ...
                                                    dp_info,...
                                                    'oil',   fimodel.oil, ...
                                                    'water', fimodel.water);
            transportModel = TransportOilWaterDPModel(G, rock, fluid, ...
                                                    rock_matrix, fluid_matrix, ....
                                                    dp_info,...
                                                    'oil',   fimodel.oil, ...
                                                    'water', fimodel.water);
                                                
            pressureModel.transfer_model_object = fimodel.transfer_model_object;
            transportModel.transfer_model_object = fimodel.transfer_model_object;
            
            pressureModel.transfer_model_name = fimodel.transfer_model_name;
            transportModel.transfer_model_name = fimodel.transfer_model_name;
       
        otherwise
            error('mrst:getSequentialModelFromFI', ...
            ['Sequential model not implemented for ''' class(fimodel), '''']);
    end
    pressureModel.operators = fimodel.operators;
    transportModel.operators = fimodel.operators;
    
    model = SequentialPressureTransportDPModel(pressureModel, transportModel, varargin{:});
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

