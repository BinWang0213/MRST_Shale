classdef DualPorosityReservoirModel < ReservoirModel
    %Base class for physical models
    %
    % SYNOPSIS:
    %   model = ReservoirModel(G, rock, fluid)
    %
    % DESCRIPTION:
    %   Extension of PhysicalModel class to accomodate reservoir-specific
    %   features such as fluid and rock as well as commonly used phases and
    %   variables.
    %
    % REQUIRED PARAMETERS:
    %   G     - Simulation grid.
    %
    %   rock  - Valid rock used for the model.
    %
    %   fluid - Fluid model used for the model.
    %
    %
    % OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
    %   See class properties.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, PhysicalModel

    properties    
        % Fluid model for the matrix
        fluid_matrix
        % Rock model for the matrix
        rock_matrix

        % Dual Porosity info - useful to pass arguments that don't change 
        % to the transfer function
        dp_info

        % Transfer Model 
        transfer_model_name
        transfer_model_object

        % Names of matrix variables
        matrixVarNames
    end

    methods
        % --------------------------------------------------------------------%
        function model = DualPorosityReservoirModel(G, varargin)
            model = model@ReservoirModel(G, [], []);

            if nargin == 1 || ischar(varargin{1})
                % We were given only grid + any keyword arguments
                doSetup = false;
            else
                assert(nargin >= 3)
                % We are being called in format
                % ReservoirModel(G, rock, fluid, ...)
                model.rock  = varargin{1};
                model.fluid = varargin{2};

                % Matrix properties
                model.rock_matrix = varargin{3};
                model.fluid_matrix = varargin{4};

                % Dual Porosity Info
                model.dp_info = varargin{5};

                % Rest of arguments should be keyword/value pairs.
                varargin = varargin(6:end);

                opt = struct('transfer_model_name', 'TransferFunction');
                [opt,unparsed] = merge_options(opt, varargin{:});

                % Transfer Function
                transfer_model_handle = str2func(opt.transfer_model_name);
                try
                    model.transfer_model_object = transfer_model_handle();
                catch
                    disp(['Could not find transfer model: ', opt.transfer_model_name ' or ' 
                          'your transfer model constructor has an input argument. In '
                          'this case you should initialize your model without the "transfer_model_name" '
                          'argument and redefine the "transfer_model" field after the model definition.'])
                    disp('Using dummy transfer function (no transfer at all !)')
                    transfer_model_handle = str2func('TransferFunction');
                    model.transfer_model_object = transfer_model_handle();
                end

                % We have been provided the means, so we will execute setup
                % phase after parsing other inputs and defaults.
                doSetup = true;
            end

            model.dpMaxRel = inf;
            model.dpMaxAbs = inf;

            model.minimumPressure = -inf;
            model.maximumPressure =  inf;

            model.dsMaxAbs = .2;

            model.nonlinearTolerance = 1e-6;
            model.inputdata = [];

            model.useCNVConvergence = false;
            model.toleranceCNV = 1e-3;
            model.toleranceMB = 1e-7;

            model.matrixVarNames = {'pwm','pom','pgm','swm','som','sgm'};

            model.extraStateOutput = false;
            model.extraWellSolOutput = true;
            model.outputFluxes = true;

            % Gravity defaults to the global variable
            model.gravity = gravity();
            % TO DO HERE
            [model, unparsed] = merge_options(model, varargin{:}); %#ok

            % Base class does not support any phases
            model.water = false;
            model.gas = false;
            model.oil = false;

            if doSetup
                if isempty(G) || isempty(model.rock)
                    warning('mrst:DualPorosityReservoirModel', ...
                        'Invalid grid/rock pair supplied. Operators have not been set up.')
                else
                    model.operators = setupOperatorsTPFA(G, model.rock, 'deck', model.inputdata);
                    operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, 'deck', model.inputdata);
                    model.operators.pv_matrix = operators_matrix.pv;
                end
            end

        end

        % --------------------------------------------------------------------%
        function state = validateState(model, state)
            % Check parent class
            state = validateState@PhysicalModel(model, state);
            active = model.getActivePhases();
            nPh = nnz(active);
            nc = model.G.cells.num;
            model.checkProperty(state, 'Pressure', [nc, 1], [1, 2]);
            if nPh > 1
                model.checkProperty(state, 'Saturation', [nc, nPh], [1, 2]);
            end
        end

        % --------------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case 'wellsol'
                    % Use colon to get all variables, since the wellsol may
                    % be empty
                    index = ':';
                    fn = 'wellSol';
                case {'swm', 'matrixwater'}
                    index = 1;
                    fn = 'swm';
                case {'som', 'matrixoil'}
                    index = 1;
                    fn = 'som';
                case {'sgm', 'matrixgas'}
                    index = 1;
                    fn = 'sgm';   
                case {'pwm', 'matrix_water_pressure'}
                    index = 1;
                    fn = 'pwm';    
                case {'pom', 'matrix_oil_pressure'}
                    index = 1;
                    fn = 'pom';
                case {'pgm', 'matrix_gas_pressure'}
                    index = 1;
                    fn = 'pgm';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@ReservoirModel(model, name);
            end
        end

        % --------------------------------------------------------------------%
        function varargout = evaluateRelPerm(model, sat, varargin)
            % Evaluate the fluid relperm. Depending on the active phases,
            % we must evaluate the right fluid relperm functions and
            % combine the results. This function calls the appropriate
            % static functions.

            opt = struct('medium', 'matrix');
            [opt,unparsed] = merge_options(opt, varargin{:});
            if(strcmp(opt.medium,'matrix'))
                fluid = model.fluid_matrix;
            else
                fluid = model.fluid;
            end
            varargin = {};

            active = model.getActivePhases();
            nph = sum(active);
            assert(nph == numel(sat), ...
            'The number of saturations must equal the number of active phases.')
            varargout = cell(1, nph);
            names = model.getPhaseNames();

            if nph > 1
                fn = ['relPerm', names];
                [varargout{:}] = model.(fn)(sat{:}, fluid, varargin{:});
            elseif nph == 1
                % Call fluid interface directly if single phase
                varargout{1} = fluid.(['kr', names])(sat{:}, varargin{:});
            end
        end
    end

    methods (Static)
        % --------------------------------------------------------------------%
        function [krW, krO, krG] = relPermWOG(sw, so, sg, f, varargin)
            % Three phase, water / oil / gas relperm.
            swcon = 0;
            if isfield(f, 'sWcon')
                swcon = f.sWcon;
            end
            swcon = min(swcon, double(sw)-1e-5);

            d  = (sg+sw-swcon);
            ww = (sw-swcon)./d;
            krW = f.krW(sw, varargin{:});

            wg = 1-ww;
            krG = f.krG(sg, varargin{:});

            krow = f.krOW(so, varargin{:});
            krog = f.krOG(so,  varargin{:});
            krO  = wg.*krog + ww.*krow;
        end

        % --------------------------------------------------------------------%
        function [krW, krO] = relPermWO(sw, so, f, varargin)
            % Two phase oil-water relperm
            krW = f.krW(sw, varargin{:});
            if isfield(f, 'krO')
                krO = f.krO(so, varargin{:});
            else
                krO = f.krOW(so, varargin{:});
            end
        end

        % --------------------------------------------------------------------%
        function [krO, krG] = relPermOG(so, sg, f, varargin)
            % Two phase oil-gas relperm.
            krG = f.krG(sg, varargin{:});
            if isfield(f, 'krO')
                krO = f.krO(so, varargin{:});
            else
                krO = f.krOG(so, varargin{:});
            end
        end

        % --------------------------------------------------------------------%
        function [krW, krG] = relPermWG(sw, sg, f, varargin)
            % Two phase water-gas relperm
            krG = f.krG(sg, varargin{:});
            krW = f.krW(sw, varargin{:});
        end

        % --------------------------------------------------------------------%
        function ds = adjustStepFromSatBounds(s, ds)
            % Ensure that cellwise increment for each phase is done with
            % the same length, in a manner that avoids saturation
            % violations.
            tmp = s + ds;

            violateUpper =     max(tmp - 1, 0);
            violateLower = abs(min(tmp    , 0));

            violate = max(violateUpper, violateLower);

            [worst, jj]= max(violate, [], 2);

            bad = worst > 0;
            if any(bad)
                w = ones(size(s, 1), 1);
                for i = 1:size(s, 2)
                    ind = bad & jj == i;
                    dworst = abs(ds(ind, i));

                    w(ind) = (dworst - worst(ind))./dworst;
                end
                ds(bad, :) = bsxfun(@times, ds(bad, :), w(bad, :));
            end
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

