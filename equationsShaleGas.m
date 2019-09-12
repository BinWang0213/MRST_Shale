function [problem, state] = equationsShaleGas(state0, state, model, dt, drivingForces, varargin)
% Generate linearized problem for the single-phase ShaleGas model-modified for
% shale gas
%
% SYNOPSIS:
%   [problem, state] = equationsShaleGas(state0, state, model, dt, drivingForces)
%
% DESCRIPTION:
%   This is the core function of the single-phase ShaleGas solver with
%   black-oil style properties. This function assembles the residual
%   equations for the conservation of ShaleGas and oil as well as required
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
%   model     - ShaleGasModel-derived class. Typically,
%               equationsShaleGas will be called from the class
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
%               for the ShaleGas pressure, as well as well equations specified
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
        [p0, wellVars0{:}] = initVariablesADI(p0, sG0, wellVars0{:});  %#ok
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
% ShaleGas props (calculated at oil pressure OK?)
bG     = f.bG(p);
rhoG   = bG.*f.rhoGS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoGf  = s.faceAvg(rhoG);
mobG   = trMult./f.muG(p);
dpG     = s.Grad(p) - rhoGf.*gdz;

% Upwind properties
upcg = (double(dpG)<=0);
vG = - s.faceUpstr(upcg, mobG).*trans.*dpG;
%Gas flow slippage and diffusion
if isfield(f,'kG_app')
   apparentK=f.kG_app(p);
   vG=s.faceUpstr(upcg, apparentK).*vG;
end
% Matrix fractuer closure, pressure-dependent matrix perm
if isfield(f,'k_gangi')
   GangiK=f.k_gangi(p);
   vG=s.faceUpstr(upcg, GangiK).*vG;
end
%Explicit fracture closure for geomechanics effect
if isfield(f,'k_hydraulicfrac')
   hf_K=f.k_hydraulicfrac(p);
   vG=s.faceUpstr(upcg, hf_K).*vG;

end
%Explicit fracture closure for geomechanics effect
if isfield(f,'k_naturalfrac')
   nf_K=f.k_naturalfrac(p);
   vG=s.faceUpstr(upcg, nf_K).*vG;
end
%non-darcy flow for hydarulic fractures
if isfield(f,'k_nondarcy')
    NondarcyK=f.k_nondarcy(p);
    [NondarcyK_face, NondarcyK]=s.splitFaceCellValue(s,upcg,NondarcyK);
    B=NondarcyK_face.*trans.*abs(dpG).*1000;
    F_ND=2.0./(1+(1+4.0.*B).^0.5);
    %max(B.val)
    %min(F_ND.val)
    vG=vG.*F_ND;
end
bGvG = s.faceUpstr(upcg, bG).*vG;


if model.outputFluxes
    state = model.storeFluxes(state,[], [], vG);
end

if model.extraStateOutput
    state = model.storebfactors(state, [], [], bG);
    state = model.storeMobilities(state, [], [],mobG);
    state = model.storeUpstreamIndices(state, [], [], upcg);
end
% EQUATIONS ---------------------------------------------------------------
names = {'gas'};
types = {'cell'};


% Single Phase equation
eqs{1} = (pv/dt).*( pvMult.*bG - pvMult0.*f.bG(p0) ) + s.Div(bGvG);

% Adsorption effect
if isfield(f,'mG_ad')
    eqs{1}=eqs{1}+ (gv/dt)./f.rhoGS.*( f.mG_ad(p) - f.mG_ad(p0) );
end

% Dummy saturation
sG = ones(model.G.cells.num, 1);
[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, {sG}, {mobG}, {rhoG}, ...
                                                                 {}, {}, ...
                                                                 drivingForces);

% well equations % I changed rhoW to bW
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap,...
        p, {mobG}, {rhoG}, {}, {}, dt, opt);
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






