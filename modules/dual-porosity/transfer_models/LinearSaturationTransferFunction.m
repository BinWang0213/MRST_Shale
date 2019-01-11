classdef LinearSaturationTransferFunction < TransferFunction
 
    properties
        beta
        snmax
        snfmin
        use_smooth_activation
    end

    methods
        function transferfunction = LinearSaturationTransferFunction(beta,snmax,snfmin,act)
            transferfunction = transferfunction@TransferFunction();
            transferfunction.nphases = 2;
            
            %% Dummy shape factor.
            if(nargin<1)
                transferfunction.beta = 1e-10;
                transferfunction.snmax = 1;
                transferfunction.snfmin = 0.1;
            else
                if(nargin<2)
                    transferfunction.beta = beta;
                    transferfunction.snmax = 1;
                    transferfunction.snfmin = 0.1;
                    transferfunction.use_smooth_activation = 0;
                else
                    if(nargin<3)
                        transferfunction.beta = beta;
                        transferfunction.snmax = snmax;
                        transferfunction.snfmin = 0.1;
                        transferfunction.use_smooth_activation = 0;
                    else
                        if(nargin<4)
                            transferfunction.beta = beta;
                            transferfunction.snmax = snmax;
                            transferfunction.snfmin = snfmin;
                            transferfunction.use_smooth_activation = 0;
                        else
                            transferfunction.beta = beta;
                            transferfunction.snmax = snmax;
                            transferfunction.snfmin = snfmin;
                            transferfunction.use_smooth_activation = act;
                        end
                    end
                end
            end
            
        end
        
        function [Talpha] = calculate_transfer(ltf,model,fracture_fields,matrix_fields)
           
            %% All calculate_transfer method should have this call. This is a "sanity check" that
            % ensures that the correct structures are being sent to calculate the transfer
            ltf.validate_fracture_matrix_structures(fracture_fields,matrix_fields);
            
            %% Minimum saturation for transfer
            satmin = ltf.snfmin;
            
            %% Pressures % water saturations
            % Matrix properties
            swm = matrix_fields.swm;
            swf = fracture_fields.swf;
            snf = 1-swf;
            snm = 1-swm;
            beta = ltf.beta;
            snmax = ltf.snmax;
            phim = model.rock_matrix.poro(:,1);
            
            %% Smooth activation function
            if(ltf.use_smooth_activation)
                kf = min(model.rock.perm,[],2);
%                 kf = model.rock.perm(:,1);
                km = model.rock_matrix.perm(:,1);
                f = (1-exp(-sqrt(kf./km).*snf))./(1-exp(-sqrt(kf./km)));
%                 f = (1-exp(-swf))./(1-exp(-1));
                satmin = 0.01;
            else
                f = 1;
            end
            
            %% Transfer
%             to = f .* beta .* phim .* double(snf > satmin) .*(snmax-snm);
            to = f .* beta .* phim .*(snmax-snm);
            tw = -to;
            
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
