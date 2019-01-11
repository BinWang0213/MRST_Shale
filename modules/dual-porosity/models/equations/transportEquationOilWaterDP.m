function [problem, state] = transportEquationOilWaterDP(state0, state, model, dt, drivingForces, varargin)
% Transport equations
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'history', [],...
             'solveForWater', false, ...
             'solveForOil', true, ...
             'iteration', -1, ...
             'stepOptions', []);  % Compatibility only

opt = merge_options(opt, varargin{:});
W = drivingForces.W;

s = model.operators;
f = model.fluid;

assert(~(opt.solveForWater && opt.solveForOil));

[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');

[p0, sW0] = model.getProps(state0, 'pressure', 'water');

% Matrix properties
[pom,swm] = model.getProps(state, 'pom','swm');
[pom0,swm0] = model.getProps(state0, 'pom','swm');

% If timestep has been split relative to pressure, linearly interpolate in
% pressure.
pFlow = p;
if isfield(state, 'timestep')
    dt_frac = dt/state.timestep;
    p = p.*dt_frac + p0.*(1-dt_frac);
	pom = pom.*dt_frac + pom0.*(1-dt_frac);
end
%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [sW,swm] = initVariablesADI(sW,swm);
    else
        assert(0, 'Backwards solver not supported for splitting');
    end
end
primaryVars = {'sW','swm'};

clear tmp

% -------------------------------------------------------------------------
sO = 1 - sW;
[krW, krO] = model.evaluateRelPerm({sW, sO});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; krO = mobMult.*krO;

% Compute transmissibility
T = s.T.*transMult;

% Gravity gradient per face
gdz = model.getGravityGradient();

% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);

% Evaluate oil properties
[vO, bO, mobO, rhoO, pO, upco, dpO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz);


% Properties for Matrix
pvMultm = pvMult;
pvMultm0 = pvMult0;

% Using capillary pressure information
pcOWm = 0;
pcOWm0 = 0;
if isfield(model.fluid_matrix, 'pcOW') && ~isempty(swm)
    pcOWm  = model.fluid_matrix.pcOW(swm);
    pcOWm0  = model.fluid_matrix.pcOW(swm0);
end

pwm = pom - pcOWm;
pwm0 = pom0 - pcOWm0;

% SMALL TO DO HERE: WE USE THE REL PERMS OF THE MATRIX TO EVALUATE THE
% EFFECTIVE PERMEABILITY
som = 1-swm;
som0 = 1-swm0;
[krWm, krOm] = model.evaluateRelPerm({swm, som});

bWm = f.bW(pwm);
bOm = f.bO(pom);

bWm0 = f.bW(pwm0);
bOm0 = f.bO(pom0);

% Transfer
vb = model.G.cells.volumes;

matrix_fields.pom = pom;
matrix_fields.swm = swm;
fracture_fields.pof = p;
fracture_fields.swf = sW;

transfer_model = model.transfer_model_object;

[Talpha] = transfer_model.calculate_transfer(model,fracture_fields,matrix_fields);

Twm = vb.*Talpha{1};
Tom = vb.*Talpha{2};

gp = s.Grad(p);
Gw = gp - dpW;
Go = gp - dpO;

if model.extraStateOutput
    state = model.storebfactors(state, bW, bO, []);
    state = model.storeMobilities(state, mobW, mobO, []);
end

if ~isempty(W)
    wflux = sum(vertcat(wellSol.flux), 2);
    perf2well = getPerforationToWellMapping(W);
    wc = vertcat(W.cells);
    
    mobWw = mobW(wc);
    mobOw = mobO(wc);
    totMobw = mobWw + mobOw;

    f_w_w = mobWw./totMobw;
    f_o_w = mobOw./totMobw;

    isInj = wflux > 0;
    compWell = vertcat(W.compi);
    compPerf = compWell(perf2well, :);

    f_w_w(isInj) = compPerf(isInj, 1);
    f_o_w(isInj) = compPerf(isInj, 2);

    bWqW = bW(wc).*f_w_w.*wflux;
    bOqO = bO(wc).*f_o_w.*wflux;

    % Store well fluxes
    wflux_O = double(bOqO);
    wflux_W = double(bWqW);
    
    for i = 1:numel(W)
        perfind = perf2well == i;
        state.wellSol(i).qOs = sum(wflux_O(perfind));
        state.wellSol(i).qWs = sum(wflux_W(perfind));
    end

end

% Get total flux from state
flux = sum(state.flux, 2);
vT = flux(model.operators.internalConn);

% Stored upstream indices
[flag_v, flag_g] = getSaturationUpwind(model.upwindType, state, {Gw, Go}, vT, s.T, {mobW, mobO}, s.faceUpstr);
flag = flag_v;

upcw  = flag(:, 1);
upco  = flag(:, 2);

upcw_g = flag_g(:, 1);
upco_g = flag_g(:, 2);

mobOf = s.faceUpstr(upco, mobO);
mobWf = s.faceUpstr(upcw, mobW);

totMob = (mobOf + mobWf);
    
mobWf_G = s.faceUpstr(upcw_g, mobW);
mobOf_G = s.faceUpstr(upco_g, mobO);
mobTf_G = mobWf_G + mobOf_G;
f_g = mobWf_G.*mobOf_G./mobTf_G;
if opt.solveForWater
    f_w = mobWf./totMob;
    bWvW   = s.faceUpstr(upcw, bW).*f_w.*vT + s.faceUpstr(upcw_g, bO).*f_g.*s.T.*(Gw - Go);

	% water fracture:
    wat_fracture = (s.pv/dt).*(pvMult.*bW.*sW       - pvMult0.*f.bW(p0).*sW0    ) + s.Div(bWvW);
	wat_fracture = wat_fracture + Twm;
    % water matrix 
	wat_matrix = (s.pv_matrix/dt).*( pvMultm.*bWm.*swm - pvMultm0.*bWm0.*swm0 );
	wat_matrix = wat_matrix - Twm;
	
	if ~isempty(W)
        wat_fracture(wc) = wat_fracture(wc) - bWqW;
    end

    eqs{1} = wat_fracture;
    oil_fracture = [];
    names = {'water_fracture','water_matrix'};
    types = {'cell','cell'};
else
    f_o = mobOf./totMob;
    bOvO   = s.faceUpstr(upco, bO).*f_o.*vT + s.faceUpstr(upco_g, bO).*f_g.*s.T.*(Go - Gw);

	% oil fracture:
    oil_fracture = (s.pv/dt).*( pvMult.*bO.*(1-sW) - pvMult0.*f.bO(p0).*(1-sW0) ) + s.Div(bOvO);
	oil_fracture = oil_fracture + Tom;
	% oil matrix
	oil_matrix = (s.pv_matrix/dt).*( pvMultm.*bOm.*som - pvMultm0.*bOm0.*som0 );
	oil_matrix = oil_matrix - Tom;
	
    if ~isempty(W)
        oil_fracture(wc) = oil_fracture(wc) - bOqO;
    end
    wat_fracture = [];
    eqs{1} = oil_fracture;
    names = {'oil_fracture','oil_matrix'};
    types = {'cell','cell'};
end

tmpEqs = {wat_fracture, oil_fracture};
tmpEqs = addFluxesFromSourcesAndBC(model, tmpEqs, ...
                                   {pFlow, pFlow},...
                                   {rhoW, rhoO},...
                                   {mobW, mobO}, ...
                                   {bW, bO},  ...
                                   {sW, sO}, ...
                                   drivingForces);
if opt.solveForWater
    eqs{1} = tmpEqs{1};
	eqs{2} = wat_matrix;
else
    eqs{1} = tmpEqs{2};
	eqs{2} = oil_matrix;
end
if ~model.useCNVConvergence
    eqs{1} = eqs{1}.*(dt./s.pv);
	eqs{2} = eqs{2}.*(dt./s.pv_matrix);
end
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
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
