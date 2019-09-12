
    clear 
    clc
    close all
    opt = struct('nkr',        1, ...
                 'shouldPlot', 1 );
%     opt = merge_options(opt, varargin{:});
    
    % Load necessary modules, etc 
    mrstModule add hfm;             % hybrid fracture module
    mrstModule add ad-core;         
    mrstModule add ad-props ad-core 
    mrstModule add mrst-gui;        % plotting routines
    mrstModule add compositional; %mrstModule add ad-blackoil;
    mrstModule add upr;

    tol=1e-5;


%     nx = 10;
%     ny = 5;
%     Lx = 400;
%     Ly = 200;

    % Set up grid and rock
    x_max = 120;
    x_min = 0;
    max_nLogGrids  = 100;
    wf = 10*milli*meter;
    deltaR  = wf/2;  % Frac width will be 2*deltaR    0.001meters = 1mm 
    factor  = 1.3;

    fracLoc = 10;

    pt = zeros(1,10);
    pt(1) = deltaR;
    fprintf('Pt[%d]= %24.16f \n',1,pt(1));
    for ilogSpacnNum = 2:max_nLogGrids 
       deltaR = deltaR * factor;
       pt(ilogSpacnNum) = pt(ilogSpacnNum-1) + deltaR;
       if(pt(ilogSpacnNum) >= ((x_max-x_min)/2.0))
          fprintf('num of log grids is: %d \n', ilogSpacnNum);
          break
       end
       fprintf('Pt[%d]= %24.16f \n',ilogSpacnNum,pt(ilogSpacnNum));
    end

    ilogSpacnNum = ilogSpacnNum - 1;
    pt = pt(1:ilogSpacnNum);

    pt_under = pt(1:24);
    pt_above = pt(1:28);
    myX1 = [0, fracLoc + [-fliplr(pt_under) ,pt_above] ];
    myX = [myX1, fliplr(80-myX1)]; 

    % G = cartGrid([nx ny],[Lx Ly]);
    % G = makeLayeredGrid(G,diff(myX));

    

%     G_matrix = tensorGrid(0:200:400, 0:100:200,myX);
    
    G_matrix = tensorGrid(-2:0.5:3, -2:0.5:3,-3:0.5:3);
    
    G_matrix = computeGeometry(G_matrix);
    
    
    if (opt.shouldPlot)
        %Plot triangular grid
        figure(1)
        plotGrid(G_matrix), view(5,45)
    end

    % epsL = 0.7;

    % rock = makeRock(G, 1*micro*darcy, 0.054, 'epsL', epsL);
    G_matrix.rock=makeRock(G_matrix,100*nano*darcy,0.063);

%     %Calculating frac cells and Setting frac properties
%     [fracIndx_X,fracIndx_Y,fracIndx_Z] = meshgrid(1:G_matrix.cartDims(1), 1:G_matrix.cartDims(2), [26,80]);
%     fraccells = sub2ind(G_matrix.cartDims, reshape(fracIndx_X,numel(fracIndx_X),1),reshape(fracIndx_Y,...
%         numel(fracIndx_X),1), reshape(fracIndx_Z,numel(fracIndx_X),1));
%     G_matrix.rock.poro(fraccells) = 0.33*3/100;
%     G_matrix.rock.perm(fraccells) = 50*darcy;





%     % Globals;                                                                        % globals stuff
%     rng(0); 
%     tol4domain = tol*1e3;
%     set1 = Field(DFN('dim',3,'n',5000,'dir',0,'ddir',0,'minl',0.3,...
%                'mu',5,'maxl',10,'bbx',[tol4domain,tol4domain,tol4domain,400-tol4domain,200-tol4domain,80-tol4domain],'dip',90,'ddip',0,...
%                'shape','l','q',4),'Poly'); 
%     set2 = Field(DFN('dim',3,'n',20000,'dir',45,'ddir',-1e9,'minl',0.4,...
%                 'mu',8,'maxl',25,'bbx',[tol4domain,tol4domain,tol4domain,400-tol4domain,200-tol4domain,80-tol4domain],'dip',45,'ddip',-1e9,...
%                 'shape','l','q',4),'Poly');  
%     set3 = Field(DFN('dim',3,'n',20000,'dir',-45,'ddir',-1e9,'minl',0.4,...
%                 'mu',8,'maxl',25,'bbx',[tol4domain,tol4domain,tol4domain,400-tol4domain,200-tol4domain,80-tol4domain],'dip',-45,'ddip',-1e9,...
%                 'shape','l','q',4),'Poly');        
%     set2 = [set2;set3];
% 
%     [set1_] = processStochFracs(set1);
%     [set2_] = processStochFracs(set2); 
%     
%     fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
%         numel(set1_),numel(set1),numel(set2_),numel(set2));
%     if (opt.shouldPlot)
% %         figure(1)
% %         Draw('ply',set1);
% %         Draw('ply',set2);view(45,30)
%         
%         figure(1)
%         Draw('ply',set1_);
%         Draw('ply',set2_);view(45,30)
%     end         
% 
%     fracSet = [set1_ ;set2_]; 
%     
%     fracplanes = struct;
%     
%     for i=1:numel(fracSet)
%         fracplanes(i).points = fracSet{i}(1:end-1,:);
%         fracplanes(i).aperture = 1*micro*meter; %1*micro*meter;
%         fracplanes(i).poro=0.5;
%         fracplanes(i).perm=100*micro*darcy;%0.01*darcy;
%     end

% # The fractures are specified by their vertices, stored in a numpy array
% f_1 = pp.Fracture(np.array([[0, 1, 2, 0], [0, 0, 1, 1], [0, 0, 1, 1]]))
% f_2 = pp.Fracture(np.array([[0.5, 0.5, 0.5, 0.5], [-1, 2, 2, -1], [-1, -1, 2, 2]]))    

    fracplanes = struct;
    fracplanes(1).points = [0 0 0;
                            1 0 0;
                            2 1 1;
                            0 1 1];
    fracplanes(1).aperture = 3*micro*meter;
    fracplanes(1).poro=0.8;
    fracplanes(1).perm=0.01*darcy;
    fracplanes(1).epsL=1.0;

    fracplanes(2).points = [0.5 -1 -1;
                            0.5  2 -1;
                            0.5  2  2;
                            0.5 -1  2];
    fracplanes(2).aperture = 3*micro*meter;
    fracplanes(2).poro=0.8;
    fracplanes(2).perm=0.01*darcy;
    fracplanes(2).epsL=1.0;

    fem = iscoplanar(fracplanes(2).points);
    if ~all(fem) 
        fprintf('The frac is not coplanar');
    end
    
    checkIfCoplanar(fracplanes)

    if (opt.shouldPlot)
        plotfracongrid(G_matrix,fracplanes); % visualize to check before pre-process
    end
    G=G_matrix;

    
%     Boi = 1.274424026342976;
%     STOIIP_NF = 6.289811*sum(G_matrix.cells.volumes .* G_matrix.rock.poro)*(1-0.23)/Boi %in stb
%     STOIIP = 6.289811*sum(G.cells.volumes .* G.rock.poro)*(1-0.23)/Boi %in stb
    % GlobTri = globalTriangulation(G_matrix);
    % [G,fracplanes] = preProcessingFractures(G, fracplanes, ...
    %                  'GlobTri', GlobTri);
    %% Process fracture(s)
    [G,fracplanes]=EDFMgrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));


    
    
    %% Fracture-Matrix NNCs
    G=fracturematrixNNC3D(G,tol);


    %%
    [G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);

    %%
    % MRST includes both natural variables and overall composition. This toggle
    % can switch between the modes.
    useNatural = true;


    % Name of problem and pressure range
    casename = 'bakken';
    pwf = 2000*psia;
    pinj = 8000*psia;


    % Set up EDFM operators
    TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);

    %% Define fluid properties

%     Pc = 170*barsa;
%     alpha = 0.5;
%     Pmax = 200*barsa;
%     m = 0.5;

    %% Define three-phase compressible flow model
    % We define a three-phase black-oil model without dissolved gas or vaporized
    % oil. This is done by first instantiating the blackoil model, and then
    % manually passing in the internal transmissibilities and the topological
    % neighborship from the embedded fracture grid.
    % gravity reset off
    % model = ThreePhaseBlackOilModel(G, [], fluid, 'disgas', false, 'vapoil', false);
    % model.operators = TPFAoperators;

    [fluid, info] = getCompositionalFluidCase(casename);
    
    eosname = 'prcorr';  %'srk','rk','prcorr'
    G1cell = cartGrid([1 1],[1 1]);
    G1cell = computeGeometry(G1cell);
    EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%     flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, 800, 10]);
%     model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', false);

    %Surface Conditions
    p_sc = 101325; %atmospheric pressure
    T_sc = 288.706;% 60 Farenheit
    [~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);


%     nkr = 1;
    flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);

    gravity reset on

    if useNatural
        model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', true);
    else
        model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', true);
    end
    model.operators = TPFAoperators;



    %% Set up initial state and schedule
    % We set up a initial state with the reference pressure and a mixture of
    % water and oil initially. We also set up a simple-time step strategy that
    % ramps up gradually towards 30 day time-steps.

    %% Set up initial state and schedule
    % We set up a initial state with the reference pressure and a mixture of
    % water and oil initially. We also set up a simple-time step strategy that
    % ramps up gradually towards 30 day time-steps.

    totTime = 10*day;
    nSteps =5;

    ncomp = fluid.getNumberOfComponents();
    s0 = [0.23, 0.70, 0.07];   %s0 = [0.23, 0.77, 0.07];

    %                                 (G, p, T, s0, z0, eos)
    state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

    wellRadius = 0.01;
    W = [];
 
    W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 1, ...
        'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
    W(1).components = info.initial;


%     switch case2run
%     case 'ProdTop' 
%         % Producer
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 26, ...
%             'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
%         W(1).components = info.initial;
%     case 'ProdTop_InjBot'
%         % Producer
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 26, ...
%             'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
%         % Injector
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
%             'comp_i', [0 0 1],'Name', 'Inj_Bot', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
%         W(1).components = info.initial;
%         W(2).components = info.injection;
%     case 'ProdBot' 
%         % Producer
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
%             'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
%         W(1).components = info.initial;
%     case 'ProdBot_InjTop'
%         % Producer
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
%             'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius); 
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 26, ...
%             'comp_i', [0 0 1],'Name', 'Inj_Bot', 'Val', pinj, 'sign', 1, 'Type', 'bhp','Radius', wellRadius);
%         W(1).components = info.initial;
%         W(2).components = info.injection;
%     otherwise
%         warning('Case Does Not Exist. Running case with Prod Only at Bottom')
%         % Producer
%         W = verticalWell(W, G.Matrix, G.Matrix.rock, 1, 1, 80, ...
%             'comp_i', [0.23, 0.76, 0.01],'Name', 'Prod_Top', 'Val', pwf, 'sign', -1, 'Type', 'bhp','Radius', wellRadius);
%         W(1).components = info.initial;
%     end
%     % plotWell(G,W);
    
    


    % bc = [];
    % bc = pside(bc, G, 'LEFT', 90.0*barsa, 'sat', [0.3, 0.3, 0.4]);
    % bc.components = info.initial;

    % state  = initResSol(G, pRef, s0);
    dt = rampupTimesteps(totTime, 2*day, nSteps);
    % dt = repmat(30*day,61,1);
    % schedule = simpleSchedule(dt, 'bc', bc);
    schedule = simpleSchedule(dt, 'W', W);

    % % Manually set the injection composition
    % [schedule.control.W.components] = deal([1     0     0     0     0     0     0     0]);
    % % Injection is pure gas
    % [schedule.control.W.compi] = deal([1, 0]);
    

    %% Simulate problem
%     fn = getPlotAfterStep(state, model, schedule);
%     [ws, states, report] = simulateScheduleAD(state, model, schedule, ...
%        'afterStepFn', getPlotAfterStep(state, model, schedule));
   
   
    [ws, states] = simulateScheduleAD(state, model, schedule, 'Verbose', true);

    %% plotting
    % %% Plot the results in the interactive viewer
    % figure(2); 
    % plotToolbar(G, states)
    % view(40,30);
    % axis tight equal;
    % 
    % figure(3);
    % plotWellSols(ws,cumsum(schedule.step.val)/86400)
%     tinSecs = cumsum(schedule.step.val);
