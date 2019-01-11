function [ state, model, schedule ] = three_phase_water_alternating_gas_setup

    %% Set up grid and rock structure
    % We define a 50x50x1 grid, spanning a 1 km by 1km domain. The porosity is
    % assigned via a Gaussian field and a synthethic permeability is computed
    % using a standard correlation.

    cartDims = [50, 50, 1];
    physDims = [1000, 1000, 1]*meter;
    G = cartGrid(cartDims, physDims);
    G = computeGeometry(G);
    gravity reset off

    rng(0);
    poro = gaussianField(cartDims, [0.05, 0.3], 3, 8);
    poro = poro(:);
    perm = poro.^3 .* (5e-5)^2 ./ (0.81 * 72 * (1 - poro).^2);
    % Build rock object
    rock = makeRock(G, perm, poro);
    
    rock_fracture.poro = 0.1*rock.poro;
    rock_fracture.perm = 100*rock.perm;

    %% Set up wells and schedule
    pv = poreVolume(G, rock);
    T = 20*year;
    irate = sum(pv)/(T*4);
    % Function handle for easily creating multiple injectors
    makeInj = @(W, name, I, J, compi) verticalWell(W, G, rock_fracture, I, J, [],...
        'Name', name, 'radius', 5*inch, 'sign', 1, 'Type', 'rate',...
        'Val', irate, 'comp_i', compi);
    W = [];
    W = makeInj(W, 'I1', 1,           1,           []);
    W = makeInj(W, 'I3', cartDims(1), cartDims(2), []);
    W = makeInj(W, 'I4', 1,           cartDims(2), []);
    W = makeInj(W, 'I2', cartDims(1), 1,           []);

    I = ceil(cartDims(1)/2);
    J = ceil(cartDims(2)/2);
    % Producer
    W = verticalWell(W, G, rock_fracture, I, J, [], 'Name', 'P1', 'radius', 5*inch, ...
        'Type', 'bhp', 'Val', 100*barsa, 'comp_i', [1, 1, 1]/3, 'Sign', -1);
    % Create two copies of the wells: The first copy is set to water injection
    % and the second copy to gas injection.
    [W_water, W_gas] = deal(W);
    for i = 1:numel(W)
        if W(i).sign < 0
            % Skip producer
            continue
        end
        W_water(i).compi = [1, 0, 0];
        W_gas(i).compi   = [0, 0, 1];
    end
    dT_target = 1*year;
    dt = rampupTimesteps(T, dT_target, 10);

    % Set up a schedule with two different controls.
    schedule = struct();
    schedule.control = [struct('W', W_water);... % Water control 1
                        struct('W', W_gas)]; % Gas control 2
    % Set timesteps
    schedule.step.val = dt;
    schedule.step.control = (cumsum(dt) > T/2) + 1;

    %% Set up fluid and simulation model
    fluid = initSimpleADIFluid('phases',    'WOG', ...
                               'rho',       [1000, 700, 250], ...
                               'n',         [2, 2, 2], ...
                               'c',         [0, 1e-3, 1e-2]/barsa, ...
                               'mu',        [1, 4, 0.25]*centi*poise ...
                               );
    
    model = ThreePhaseBlackOilDPModel(G, rock_fracture, fluid,rock, fluid, [],'disgas', false, 'vapoil', false);
    % The shape factor and transfer function
    fracture_spacing = repmat([1,1,1],G.cells.num,1);
    shape_factor_name = 'KazemiShapeFactor';
    model.transfer_model_object = KazemiOilWaterGasTransferFunction(shape_factor_name,fracture_spacing);
    
    % Set up initial reservoir at 100 bar pressure and completely oil filled.
    state = initResSol(G, 100*barsa, [0, 1, 0]);
    state.pom = state.pressure;
    state.swm = state.s(:,1);
    state.som = state.s(:,2);
    state.sgm = state.s(:,3);
end

