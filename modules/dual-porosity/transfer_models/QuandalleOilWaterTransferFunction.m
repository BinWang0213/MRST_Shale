classdef QuandalleOilWaterTransferFunction < TransferFunction
 
    properties
        block_dimensions
    end
    
    methods
        
        function transferfunction = QuandalleOilWaterTransferFunction(block_dimensions)
            
            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 2;
            
            %% Dummy shape factor.
            if(nargin>0)
                transferfunction.block_dimensions = block_dimensions;
            else
                transferfunction.block_dimensions = ones(1,3);
            end
            
        end
        
        function [Talpha] = calculate_transfer(ktf,model,fracture_fields,matrix_fields)

            %% All calculate_transfer method should have this call. This is a "sanity check" that
            % ensures that the correct structures are being sent to calculate the transfer
            ktf.validate_fracture_matrix_structures(fracture_fields,matrix_fields);                                         
                                                                                  
            %% The varibles
            pom = matrix_fields.pom;
            swm = matrix_fields.swm;
            pof = fracture_fields.pof;
            swf = fracture_fields.swf;
            
            %% Oil Saturations
            %Fracture Oil saturation
            sof  = 1 - swf;
            
            %Matrix Oil saturation
            som = 1 - swm;
            
            %% Pressures 
            pcOW = 0;
            pcOWm = 0;
            
            if isfield(model.fluid, 'pcOW') && ~isempty(swf)
                pcOW  = model.fluid.pcOW(swf);
            end

            if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
                pcOWm  = model.fluid_matrix.pcOW(swm);
            end

            pwm = pom - pcOWm;
            pwf = pof - pcOW;
            
             %% Evaulate Rel Perms
            %Rel perms for the transfer
            [krwf, krnf] = model.evaluateRelPerm({swf, sof});
            [krwm, krnm] = model.evaluateRelPerm({swm, som});
            
            %% Additional Properties
            km = model.rock_matrix.perm(:,1);
            muw = model.fluid_matrix.muW(pwm);
            mun = model.fluid_matrix.muO(pom);
%             bWm = model.fluid.bW(pwm);
%             bOm = model.fluid.bO(pom);
%             bW = model.fluid.bW(pW);
%             bO = model.fluid.bO(pO);
           
            %% rho star
%             rhostar = fracture_fields.rhostar;
            rhow = model.fluid.rhoWS;
            rhoo = model.fluid.rhoOS;
            rhostar = swf*rhow+(1-swf)*rhoo;
            
            %% Shape factors on each direction
            a = ktf.block_dimensions(:,1);
            b = ktf.block_dimensions(:,2);
            c = ktf.block_dimensions(:,3);
            
            if(size(a)==[1,1])
                a = a*ones(size(model.rock_matrix.poro));
            end
            
            if(size(b)==[1,1])
                b = b*ones(size(model.rock_matrix.poro));
            end
            
            if(size(c)==[1,1])
                
                c = c*ones(size(model.rock_matrix.poro));
            end
            
            %% Gravity
            g = 9.81;
            
            %% Potentials
            % Wetting phase
            w_m_psi_z_m = pwm + g*rhow.*c/2;
            w_m_psi_z_p = pwm - g*rhow.*c/2;
            w_m_psi_x_m = pwm;
            w_m_psi_x_p = pwm;
            w_m_psi_y_m = pwm;
            w_m_psi_y_p = pwm;

            w_f_psi_z_m = pwf + g*rhostar.*c/2;
            w_f_psi_z_p = pwf - g*rhostar.*c/2;
            w_f_psi_x_m = pwf;
            w_f_psi_x_p = pwf;
            w_f_psi_y_m = pwf;
            w_f_psi_y_p = pwf;

            % Non-wetting phase
            n_m_psi_z_m = pom + g*rhoo.*c/2;
            n_m_psi_z_p = pom - g*rhoo.*c/2;
            n_m_psi_x_m = pom;
            n_m_psi_x_p = pom;
            n_m_psi_y_m = pom;
            n_m_psi_y_p = pom;

            n_f_psi_z_m = pof + g*rhostar.*c/2;
            n_f_psi_z_p = pof - g*rhostar.*c/2;
            n_f_psi_x_m = pof;
            n_f_psi_x_p = pof;
            n_f_psi_y_m = pof;
            n_f_psi_y_p = pof;

            % Potential Differences
            w_psi_diff_z_m = w_f_psi_z_m - w_m_psi_z_m;
            w_psi_diff_z_p = w_f_psi_z_p - w_m_psi_z_p;
            w_psi_diff_x_m = w_f_psi_x_m - w_m_psi_x_m;
            w_psi_diff_x_p = w_f_psi_x_p - w_m_psi_x_p;
            w_psi_diff_y_m = w_f_psi_y_m - w_m_psi_y_m;
            w_psi_diff_y_p = w_f_psi_y_p - w_m_psi_y_p;

            n_psi_diff_z_m = n_f_psi_z_m - n_m_psi_z_m;
            n_psi_diff_z_p = n_f_psi_z_p - n_m_psi_z_p;
            n_psi_diff_x_m = n_f_psi_x_m - n_m_psi_x_m;
            n_psi_diff_x_p = n_f_psi_x_p - n_m_psi_x_p;
            n_psi_diff_y_m = n_f_psi_y_m - n_m_psi_y_m;
            n_psi_diff_y_p = n_f_psi_y_p - n_m_psi_y_p;

            % Shape factors
            sigma_x = 2./a.^2;
            sigma_y = 2./b.^2;
            sigma_z = 2./c.^2;

            %% Mobilities
            %% This flags equals 1 for each cell if flow is coming from 
            % the fractures and zero otherwise. 
            flow_direction_w_xp = double(double(w_psi_diff_x_p)>=0);
            flow_direction_o_xp = double(double(n_psi_diff_x_p)>=0);

            flow_direction_w_xm = double(double(w_psi_diff_x_m)>=0);
            flow_direction_o_xm = double(double(n_psi_diff_x_m)>=0);

            flow_direction_w_yp = double(double(w_psi_diff_y_p)>=0);
            flow_direction_o_yp = double(double(n_psi_diff_y_p)>=0);

            flow_direction_w_ym = double(double(w_psi_diff_y_m)>=0);
            flow_direction_o_ym = double(double(n_psi_diff_y_m)>=0);

            flow_direction_w_zp = double(double(w_psi_diff_z_p)>=0);
            flow_direction_o_zp = double(double(n_psi_diff_z_p)>=0);

            flow_direction_w_zm = double(double(w_psi_diff_z_m)>=0);
            flow_direction_o_zm = double(double(n_psi_diff_z_m)>=0);

            krwxpt = krwf.*flow_direction_w_xp + krwm.*(~flow_direction_w_xp);
            krnxpt = krnf.*flow_direction_o_xp + krnm.*(~flow_direction_o_xp);

            krwxmt = krwf.*flow_direction_w_xm + krwm.*(~flow_direction_w_xm);
            krnxmt = krnf.*flow_direction_o_xm + krnm.*(~flow_direction_o_xm);

            krwypt = krwf.*flow_direction_w_yp + krwm.*(~flow_direction_w_yp);
            krnypt = krnf.*flow_direction_o_yp + krnm.*(~flow_direction_o_yp);

            krwymt = krwf.*flow_direction_w_ym + krwm.*(~flow_direction_w_ym);
            krnymt = krnf.*flow_direction_o_ym + krnm.*(~flow_direction_o_ym);

            krwzpt = krwf.*flow_direction_w_zp + krwm.*(~flow_direction_w_zp);
            krnzpt = krnf.*flow_direction_o_zp + krnm.*(~flow_direction_o_zp);

            krwzmt = krwf.*flow_direction_w_zm + krwm.*(~flow_direction_w_zm);
            krnzmt = krnf.*flow_direction_o_zm + krnm.*(~flow_direction_o_zm);

            mobwxp = krwxpt./muw;
            mobnxp = krnxpt./mun;

            mobwxm = krwxmt./muw;
            mobnxm = krnxmt./mun;

            mobwyp = krwypt./muw;
            mobnyp = krnypt./mun;

            mobwym = krwymt./muw;
            mobnym = krnymt./mun;

            mobwzp = krwzpt./muw;
            mobnzp = krnzpt./mun;

            mobwzm = krwzmt./muw;
            mobnzm = krnzmt./mun;
    
            %% Transfer
            Tw = sigma_x .* km .* mobwxm .* w_psi_diff_x_m + ...
                 sigma_x .* km .* mobwxp .* w_psi_diff_x_p + ...
                 sigma_y .* km .* mobwym .* w_psi_diff_y_m + ...
                 sigma_y .* km .* mobwyp .* w_psi_diff_y_p + ...
                 sigma_z .* km .* mobwzm .* w_psi_diff_z_m + ...
                 sigma_z .* km .* mobwzp .* w_psi_diff_z_p;

            To = sigma_x .* km .* mobnxm .* n_psi_diff_x_m + ...
                 sigma_x .* km .* mobnxp .* n_psi_diff_x_p + ...
                 sigma_y .* km .* mobnym .* n_psi_diff_y_m + ...
                 sigma_y .* km .* mobnyp .* n_psi_diff_y_p + ...
                 sigma_z .* km .* mobnzm .* n_psi_diff_z_m + ...
                 sigma_z .* km .* mobnzp .* n_psi_diff_z_p;
             
            Talpha{1} = Tw;
            Talpha{2} = To;
            
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
