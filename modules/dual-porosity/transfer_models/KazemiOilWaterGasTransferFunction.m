classdef KazemiOilWaterGasTransferFunction < TransferFunction
 
    properties
        shape_factor_object
        shape_factor_name
    end
    
    methods
        function transferfunction = KazemiOilWaterGasTransferFunction(shape_factor_name,fracture_spacing)
            
            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 3;
            
            if (nargin<1)
                %% No information about shape factor is provided so we set it to 1
                transferfunction.shape_factor_name = 'ConstantShapeFactor';
                fracture_spacing = [1,1,1];
                shape_factor_value = 1;
                shape_factor_handle = str2func(transferfunction.shape_factor_name);
                transferfunction.shape_factor_object = shape_factor_handle(fracture_spacing,shape_factor_value);
            else
                if (nargin<2)
                    msg = ['ERROR: If you provide the shape factor name, make sure to',...
                          'provide also the fracture spacing, [lx,ly,lz]'];
                    error(msg);
                else
                    transferfunction.shape_factor_name = shape_factor_name;
                    shape_factor_handle = str2func(transferfunction.shape_factor_name);
                    transferfunction.shape_factor_object = shape_factor_handle(fracture_spacing);
                end
            end
            
        end
        
        function [Talpha] = calculate_transfer(ktf,model,fracture_fields,matrix_fields)

            %% All calculate_transfer method should have this call. This is a "sanity check" that
            % ensures that the correct structures are being sent to calculate the transfer
            ktf.validate_fracture_matrix_structures(fracture_fields,matrix_fields);                                         
                                                                                  
            %% The varibles
            pom = matrix_fields.pom;
            swm = matrix_fields.swm;
            sgm = matrix_fields.sgm;
            
            pO = fracture_fields.pof;
            sW = fracture_fields.swf;
            sG = fracture_fields.sgf;
            
            %% Reconstructing the other phase pressures
            
            %% Pressures 
            pcOW = 0;
            pcOWm = 0;
            pcOG = 0;
            pcOGm = 0;
            
            if isfield(model.fluid, 'pcOW') && ~isempty(sW)
                pcOW  = model.fluid.pcOW(sW);
                pcOG  = model.fluid.pcOG(sG);
            end

            if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
                pcOWm  = model.fluid_matrix.pcOW(swm);
                pcOGm  = model.fluid_matrix.pcOG(sgm);
            end
            
            pwm = pom - pcOWm;
            pgm = pom + pcOGm;
            
            pW = pO - pcOW;
            pG = pO + pcOG;
            
            %% Oil Saturations
            %Fracture Oil saturation
            sO  = 1 - sW - sG;
            
            %Matrix Oil saturation
            som = 1 - swm - sgm;

            %% Evaulate Rel Perms
            %Rel perms for the transfer
            [krW, krO, krG] = model.evaluateRelPerm({sW, sO, sG});
            [krWm, krOm, krGm] = model.evaluateRelPerm({swm, som, sgm});
            
            %% This flags equals 1 for each cell if flow is coming from 
            % the fractures and zero otherwise. 
            dpw = (double(pwm-pW)<=0);
            dpo = (double(pom-pO)<=0);
            dpg = (double(pgm-pG)<=0);
            
            krwt = krW.*dpw + krWm.*(~dpw);
            krot = krO.*dpo + krOm.*(~dpo);
            krgt = krG.*dpg + krGm.*(~dpg);
            
            %% Additional Properties
            km = model.rock_matrix.perm(:,1);
            muwm = model.fluid_matrix.muW(pwm);
            muom = model.fluid_matrix.muO(pom);
            mugm = model.fluid_matrix.muG(pgm);
            
            bWm = model.fluid.bW(pwm);
            bOm = model.fluid.bO(pom);
            bGm = model.fluid.bG(pgm);
            
			
			%% Shape Factor
            [sigma]=ktf.shape_factor_object.calculate_shape_factor(model);
           
            %% Compute Transfer
            %(units 1/T)
			% Note: shape factors include permeability
            tw=(sigma.*bWm.*krwt./muwm).*(pW-pwm); 
            to=(sigma.*bOm.*krot./muom).*(pO-pom);
            tg=(sigma.*bGm.*krgt./mugm).*(pG-pgm);
            
            %% Note that we return a 3x1 Transfer since our model is 3ph
            Talpha{1} = tw;
            Talpha{2} = to;
            Talpha{3} = tg;
        end
        
        function [] = validate_fracture_matrix_structures(ktf,fracture_fields,matrix_fields)
            %% We use the superclass to validate the structures of matrix/fracture variables                                          
            validate_fracture_matrix_structures@TransferFunction(ktf,fracture_fields,matrix_fields);
        end
        
    end
    
    
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
