
    clear 
    clc
    close all
    addpath('./utils')
    addpath('../ADFNE15')
    addpath('../mrst-2018b')
    
    startup
    
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

%     input parameters
    [xmin,ymin,zmin,xmax,ymax,zmax] = deal(0,0,0,120,500,50);
    tol=1e-5;
    totTime = 20*year;
    nSteps =25;
    wellRadius = 0.01;
        

    % Set up grid and rock
    x_max = 120;
    x_min = 0;
    max_nLogGrids  = 100;
    deltaR         = 0.05;  % Frac width will be 2*deltaR    0.001meters = 1mm 
    factor         = 1.3;

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
    myX = 60 + [-fliplr(pt) ,pt];
    %increase from 5 to 20
    temp = 320+cumsum(fliplr([6.3956  8.3643   10.9236   14.2507   18.5759]));

    myY = [ymin:20:320 temp 385:5:ymax];
    
    G_matrix = tensorGrid(myY,[xmin myX xmax],[zmin zmax]);
    G_matrix = computeGeometry(G_matrix);
    
    
    if (opt.shouldPlot)
        %Plot triangular grid
        figure(1)
        plotGrid(G_matrix), view(5,45)
    end

    G_matrix.rock = makeRock(G_matrix, 1.0e-19, 0.04);  %100*nano*darcy
    %set fracture properties
    fracIndx_X = 28:G_matrix.cartDims(2);
    numFracCells = length(fracIndx_X);
    fracIndx_Y = ones(numFracCells,1)*ceil(G_matrix.cartDims(1)/2);
    fraccells = sub2ind(G_matrix.cartDims, fracIndx_Y, fracIndx_X');
    G_matrix.rock.poro(fraccells) = 0.33*3/100;
    G_matrix.rock.perm(fraccells) = 5e-11; %50*darcy;

    rng(0); 
    tol4domain = tol*1e3;
    
%   Microfractures parallel to bedding planes
    set1 = Field(DFN('dim',3,'n',8000,'dir',0,'ddir',0,'minl',0.3,...
               'mu',5,'maxl',7,'bbx',[tol4domain,tol4domain,tol4domain,ymax-tol4domain,xmax-tol4domain,zmax-tol4domain],'dip',90,'ddip',0,...
               'shape','l','q',4),'Poly'); 
    set2 = Field(DFN('dim',3,'n',8000,'dir',45,'ddir',-1e9,'minl',0.4,...
                'mu',8,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,ymax-tol4domain,xmax-tol4domain,zmax-tol4domain],'dip',45,'ddip',-1e9,...
                'shape','l','q',4),'Poly');  
    set3 = Field(DFN('dim',3,'n',8000,'dir',-45,'ddir',-1e9,'minl',0.4,...
                'mu',8,'maxl',15,'bbx',[tol4domain,tol4domain,tol4domain,ymax-tol4domain,xmax-tol4domain,zmax-tol4domain],'dip',-45,'ddip',-1e9,...
                'shape','l','q',4),'Poly');
            
    set2 = [set2;set3];

    [set1_,nonPlanarSets1,fracArea1] = processStochFracs(set1);
    [set2_,nonPlanarSets2,fracArea2] = processStochFracs(set2); 
%     fracArea = [fracArea1;fracArea2];
    
    fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
        numel(set1_),numel(set1),numel(set2_),numel(set2));
    
    fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
        numel(nonPlanarSets1),numel(nonPlanarSets2));

    fracSet = [set1_ ;set2_]; 
    
    if (opt.shouldPlot)    
        figure(2)
%         Draw('ply',set1_); view(5,45)
%         Draw('ply',set2_); view(5,45)
        Draw('ply',fracSet); view(5,45)        
    end         
    

    
    fracplanes = struct;
    
    for i=1:numel(fracSet)
        fracplanes(i).points = fracSet{i}(1:end-1,:);
        fracplanes(i).aperture = 1*micro*meter; %1*micro*meter;
        fracplanes(i).poro=0.5;
        fracplanes(i).perm=100*micro*darcy;%0.01*darcy;
    end
      
%     checkIfCoplanar(fracplanes)

    if (opt.shouldPlot)
        plotfracongrid(G_matrix,fracplanes); % visualize to check before pre-process
    end
    G=G_matrix;

    
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
    casename = 'onlymethane';
    pwf = 500*psia;
    p_init = 5000*psia;
    
    % Set up EDFM operators
    TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);

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

    %Surface Conditions
    p_sc = 101325; %atmospheric pressure
    T_sc = 288.706;% 60 Farenheit
    [~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

    flowfluid = initSimpleADIFluid('n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]);    % flowfluid.KGangiFn = @(p) power((1-power(((Pc - alpha.*p)./Pmax),m)),3);

    gravity reset on

    if useNatural
        model = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', false);
    else
        model = OverallCompositionCompositionalModel(G, [], flowfluid, fluid, 'water', false);
    end
    model.operators = TPFAoperators;


    %% Set up initial state and schedule
    % We set up a initial state with the reference pressure and a mixture of
    % water and oil initially. We also set up a simple-time step strategy that
    % ramps up gradually towards 30 day time-steps.

    ncomp = fluid.getNumberOfComponents();
    s0 = [1];  

    %                                 (G, p, T, s0, z0, eos)
    state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

    W = [];
    
    % Producer
    ind = [G.cartDims(1) ceil(G.cartDims(2)/2)];
    prodcell = sub2ind(G.cartDims, ind(1), ind(2));
    W = verticalWell(W, G.Matrix, G.Matrix.rock, prodcell, [],'Radius', 0.1, ...
        'comp_i', [1], 'Name', 'Prod', 'Val', pwf, 'sign', -1, 'Type', 'bhp');
    
    W(1).components = info.initial;

    dt =rampupTimesteps(totTime, 2*year, nSteps);
    schedule = simpleSchedule(dt, 'W', W);

    % % Manually set the injection composition
    % [schedule.control.W.components] = deal([1     0     0     0     0     0     0     0]);
    % % Injection is pure gas
    % [schedule.control.W.compi] = deal([1, 0]);
    

    %% Simulate problem
%     fn = getPlotAfterStep(state, model, schedule);
%     [ws, states, report] = simulateScheduleAD(state, model, schedule, ...
%        'afterStepFn', getPlotAfterStep(state, model, schedule));
   
   
      [ws, states] = simulateScheduleAD(state, model, schedule, 'Verbose', false);

    %% plotting
    if  opt.shouldPlot
        figure(3)
        plotToolbar(G, states);
        axis equal tight off
        daspect([1 1 0.2])
        view(85, 20);
        
        plotWellSols(ws, cumsum(schedule.step.val), 'field','qGs','linestyles',{'r'})
    end
    


