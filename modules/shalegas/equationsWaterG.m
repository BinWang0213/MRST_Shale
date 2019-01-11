function [problem, state] = equationsWaterG(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the single-phase water model-modified for
% shale gas
%
% SYNOPSIS:
%   [problem, state] = equationsWater(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the single-phase water solver with
%   black-oil style properties. This function assembles the residual
%   equations for the conservation of water and oil as well as required
%   well equations. By default, Jacobians are also provided by the use of
%   automatic differentiation.
%
% REQUIRED PARAMETERS:
%   state0    - Reservoir state at the previous timestep. Assumed to have
%               physically reasonable values.
%
%   state     - State at the current nonlinear iteration. The values do not
%               need to be physically reasonable.
%
%   model     - WaterModel-derived class. Typically,
%               equationsWater will be called from the class
%               getEquations member function.
%
%   dt        - Scalar timestep in seconds.
%
%   drivingForces - Struct with fields:
%                   * W for wells. Can be empty for no wells.
%                   * bc for boundary conditions. Can be empty for no bc.
%                   * src for source terms. Can be empty for no sources.
%
% OPTIONAL PARAMETERS:
%   'Verbose'    -  Extra output if requested.
%
%   'reverseMode'- Boolean indicating if we are in reverse mode, i.e.
%                  solving the adjoint equations. Defaults to false.
%
%   'resOnly'    - Only assemble residual equations, do not assemble the
%                  Jacobians. Can save some assembly time if only the
%                  values are required.
%
%   'iterations' - Nonlinear iteration number. Special logic happens in the
%                  wells if it is the first iteration.
% RETURNS:
%   problem - LinearizedProblemAD class instance, containing the equation
%               for the water pressure, as well as well equations specified
%               by the WellModel class.
%
%   state   - Updated state. Primarily returned to handle changing well
%             controls from the well model.
%
% SEE ALSO:
%   equationsBlackOil, ThreePhaseBlackOilModel

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);  % Compatibility only

opt = merge_options(opt, varargin{:});

W = drivingForces.W;

s = model.operators;
G = model.G;
f = model.fluid;
rock=model.rock;

[p, wellSol] = model.getProps(state, 'pressure', 'wellsol');

[p0, wellSol0] = model.getProps(state0, 'pressure', 'wellSol');
[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

%Initialization of independent variables ----------------------------------

if ~opt.resOnly,
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode,
        [p, wellVars{:}] = initVariablesADI(p, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, wellVars0{:}] = initVariablesADI(p0, sW0, wellVars0{:});  %#ok
    end
end
primaryVars = {'pressure', wellVarNames{:}};
gdz   = s.Grad(G.cells.centroids) * model.getGravityVector()';
%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end

transMult=1;
if isfield(f, 'transMult')
   transMult=f.transMult(p); 
end

trans=s.T.*transMult;
% -------------------------------------------------------------------------
%Rock propeties
pv = (rock.poro) .* G.cells.volumes;
gv = (1.-rock.poro) .* G.cells.volumes;
% water props (calculated at oil pressure OK?)
bW     = f.bW(p);
rhoW   = bW.*f.rhoWS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoWf  = s.faceAvg(rhoW);
mobW   = trMult./f.muW(p);
dpW     = s.Grad(p) - rhoWf.*gdz;

% Upwind properties
upcw = (double(dpW)<=0);
vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
%Gas flow slippage and diffusion
if isfield(f,'kG_app')
   apparentK=f.kG_app(p);
   vW=s.faceUpstr(upcw, apparentK).*vW;
end
% Matrix fractuer closure, pressure-dependent matrix perm
if isfield(f,'k_gangi')
   GangiK=f.k_gangi(p);
   vW=s.faceUpstr(upcw, GangiK).*vW;
end
%Explicit fracture closure for geomechanics effect
if isfield(f,'k_hydraulicfrac')
   hf_K=f.k_hydraulicfrac(p);
   vW=s.faceUpstr(upcw, hf_K).*vW;

end
%Explicit fracture closure for geomechanics effect
if isfield(f,'k_naturalfrac')
   nf_K=f.k_naturalfrac(p);
   vW=s.faceUpstr(upcw, nf_K).*vW;
end
%non-darcy flow for hydarulic fractures
if isfield(f,'k_nondarcy')
    NondarcyK=f.k_nondarcy(p);
    [NondarcyK_face, NondarcyK]=s.splitFaceCellValue(s,upcw,NondarcyK);
    B=NondarcyK_face.*trans.*abs(dpW).*1000;
    F_ND=2.0./(1+(1+4.0.*B).^0.5);
    %max(B.val)
    %min(F_ND.val)
    vW=vW.*F_ND;
end
bWvW = s.faceUpstr(upcw, bW).*vW;


if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, [], []);
    state = model.storeMobilities(state, mobW, [], []);
    state = model.storeUpstreamIndices(state, upcw, [], []);
end
% EQUATIONS ---------------------------------------------------------------
names = {'water'};
types = {'cell'};


% Single Phase equation
eqs{1} = (pv/dt).*( pvMult.*bW - pvMult0.*f.bW(p0) ) + s.Div(bWvW);

% Adsorption effect
if isfield(f,'mG_ad')
    eqs{1}=eqs{1}+ (gv/dt)./f.rhoGS.*( f.mG_ad(p) - f.mG_ad(p0) );
end

% Dummy saturation
sW = ones(model.G.cells.num, 1);
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, {sW}, {mobW}, {rhoW}, ...
                                                                 {}, {}, ...
                                                                 drivingForces);

% well equations % I changed rhoW to bW
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap,...
        p, {mobW}, {rhoW}, {}, {}, dt, opt);
problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end
%--------------------------------------------------------------------------

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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






