classdef GilmanKazemiOilWaterTransferFunction < TransferFunction
    % Gilman Kazemi 1988 2 Phase Transfer Function
    
    properties
        shape_factor_object
        shape_factor_name
    end
    
    methods
        function transferfunction = GilmanKazemiOilWaterTransferFunction(shape_factor_name,block_dimension)
            
            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 2;
            
            if (nargin<1)
                %% No information about shape factor is provided so we set it to 1
                transferfunction.shape_factor_name = 'ConstantShapeFactor';
                block_dimension = [1,1,1];
                shape_factor_value = 1;
                shape_factor_handle = str2func(transferfunction.shape_factor_name);
                transferfunction.shape_factor_object = shape_factor_handle(block_dimension,shape_factor_value);
            else
                if (nargin<2)
                    msg = ['ERROR: If you provide the shape factor name, make sure to',...
                          'provide also the fracture spacing, [lx,ly,lz]'];
                    error(msg);
                else
                    transferfunction.shape_factor_name = shape_factor_name;
                    shape_factor_handle = str2func(transferfunction.shape_factor_name);
                    transferfunction.shape_factor_object = shape_factor_handle(block_dimension);
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
            pO = fracture_fields.pof;
            sW = fracture_fields.swf;
            
            %% Oil Saturations
            %Fracture Oil saturation
            sO  = 1 - sW;
            
            %Matrix Oil saturation
            som = 1-swm;
            
            %% Pressures 
            pcOW = 0;
            pcOWm = 0;
            
            if isfield(model.fluid, 'pcOW') && ~isempty(sW)
                pcOW  = model.fluid.pcOW(sW);
            end

            if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
                pcOWm  = model.fluid_matrix.pcOW(swm);
            end

            pwm = pom - pcOWm;
            pW = pO - pcOW;
            
             %% Evaulate Rel Perms
            %Rel perms for the transfer
            [krW, krO] = model.evaluateRelPerm({sW, sO});
            [krWm, krOm] = model.evaluateRelPerm({swm, som});
            
            %% This flags equals 1 for each cell if flow is coming from 
            % the fractures and zero otherwise. 
            dpw = (double(pwm-pW)<=0);
            dpo = (double(pom-pO)<=0);
            
            krwt = krW.*dpw + krWm.*(~dpw);
            krot = krO.*dpo + krOm.*(~dpo);

            %% Additional Properties
            km = model.rock_matrix.perm(:,1);
            muwm = model.fluid_matrix.muW(pwm);
            muom = model.fluid_matrix.muO(pom);
            bWm = model.fluid.bW(pwm);
            bOm = model.fluid.bO(pom);
            
            %matrix fluid densities
            rhoWm = model.fluid_matrix.rhoWS;
            rhoOm= model.fluid_matrix.rhoOS;
                       
            %fracture fluid densities
            rhoWf = model.fluid.rhoWS;
            rhoOf= model.fluid.rhoOS;
            
            %Gravity
            g = 9.80665;
            
            %% Residual/connate saturations
            
            %these are set as default 0 need to be implemented properly
            %connate water
            swcm=0; 
            swcf=0;

            %residual oil in presence of water
            sorwm=0;
            sorwf=0;
            
            %% Shape Factor
            %lx,ly,lz are the inverse of the fracture spacing given
            % by the user i.e. the matrix block dimensions
            [sigma]=ktf.shape_factor_object.calculate_shape_factor(model);
           
           %% Fluid heights
           
           %vertical matrix block height
            lz=ktf.shape_factor_object.block_dimension(3);
            
            hof=(sO-sorwf)./(1-swcf-sorwf)*lz;
            hwf=(sW-swcf)./(1-swcf-sorwf)*lz;
            hom=(som-sorwm)./(1-sorwm-swcm)*lz;
            hwm=(swm-swcm)./(1-sorwm-swcm)*lz;
            
            %fracture and matrix block centroid depths
            Df=model.G.cells.centroids(:,3);
            Dm=model.G.cells.centroids(:,3);
            %% Compute Transfer
            %(units 1/T)
            % Note: shape factors include permeability
            
            to=(sigma.*bOm.*krot./muom).*(pO-pom-rhoOf.*g.*(hof-Df)-rhoOm.*g.*(hom-Dm));
            tw=(sigma.*bWm.*krwt./muwm).*(pW-pwm-rhoWf.*g.*(hwf-Df)-rhoWm.*g.*(hwm-Dm));
%             
%             double(to)
%             double(tw)
            
            %% Note that we return a 2x1 Transfer since our model is 2ph
            Talpha{1} = tw;
            Talpha{2} = to;
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
