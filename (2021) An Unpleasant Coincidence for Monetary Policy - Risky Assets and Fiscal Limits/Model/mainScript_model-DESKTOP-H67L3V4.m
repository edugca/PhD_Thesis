%% PAPER: An Unpleasant Coincidence for Monetary Policy: Risky Assets and Fiscal Limits
%
% Author: Eduardo G. C. Amaral
% Version date: June 4th 2020
%
% This routine runs the model

%% Housekeeping
clear;
clc;
close all;

% Final or Temporary
saveToFinalFolder = false;

% Paper
paperNumber = 2;

% Paths of the model
if paperNumber == 2
    pathMain    = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) An Unpleasant Coincidence for Monetary Policy - Risky Assets and Fiscal Limits', 'C:\');
    pathData    = pub_Path('/Users/Eduardo/OneDrive/MATLAB/My Library/Database/Data/Brazil', 'C:\');
end
pathTables  = [pathMain filesep 'Tables'];
pathImages  = [pathMain filesep 'Images'];
pathSaved   = [pathMain filesep 'Saved'];

% Paths of results
if saveToFinalFolder
    if paperNumber == 2
        pathImages = edu_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 2/Images', 'C:\');
        pathTables = edu_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 2/Tables', 'C:\');
    end
end

% Path of RISE
pathRISE = edu_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/RISE_toolbox-master', 'C:\');
addpath(pathRISE);

% Graph config
custom_LineWidth = {1, 1, 1, 1, 1, 1, 1};
custom_LineStyle = {'-', '--', ':', '-.', '-', '--', ':'};
custom_LineColor = {'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red'};
custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};

% Add paths of auxiliary routines
addpath(genpath('Auxiliary routines')); % Routines also udeful for calculating the fiscal limits
addpath(genpath('Auxiliary routines 2')); % Routines useful for running the regime-switching model

% Rise start-up
rise_startup();

%% Configure the model

% Load parameters from estimation
loadEstimParameters = 'All'; % 'None', 'OnlyShocks', 'All'

% Simplify policy rule
conf_simplifyPolicyRule = false;

% Append welfare equations
conf_appendWelfareEquations = true;

% Sticky prices
conf_stickyPrices = true;

% Monopolistic competition
conf_monopolisticCompetition = true;

% Reestimate logistic approximation of the fiscal limit
reestimateFiscalLimitWithApproximation = true;

% Approximate default probability
conf_approximateDefaultProb = 'logistic'; % 'linear', 'logistic'

% Model file
%ssFileAddress   = 'paper2_model_steadyState.rs';

% Calibration name
if paperNumber == 2
    conf_calibrationName = 'globalsolutionModel'; % 'textbook', 'Brazil', 'globalsolutionModel'
elseif paperNumber == 3
    conf_calibrationName = 'globalsolutionModel'; % 'textbook', 'Brazil', 'globalsolutionModel'
end

% Occasionally binding constraint (should rename this)
conf_occBinConstraint = 'FiscalLimit'; % 'None', 'FiscalLimit'
conf_regimeSwitching = 'Endogenous'; % 'Exogenous', 'Endogenous'

% Effective lower bound
conf_effectiveLowerBound = false;

% Has government?
conf_hasGovernment = true; % true, false

% Government bonds are risk-free?
conf_govBondsAreRiskFree = false; % true, false

% Government expenses
conf_govExpenses = 'Actual'; % 'Actual', 'High', 'Zero'

% Debt level
conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
conf_defPolTarget = 0;
conf_debtLevelGuess = 2.5;

% Government accumulates debt
conf_govAccumDebt = true; % true, false

% Tax rule depends on B_t/B or B_t/Y_t * Y/B
conf_taxRule = 'Relative'; % 'Absolute', 'Relative'

% Capital status
conf_kVariable = 'Fixed'; % 'Fixed', 'Exogenous'

% Productivity process type
conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'

% DefR process type
conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'

% Policy rule
conf_policyRule = 'fixedIntercept_inflation'; % 'fixedIntercept_inflation', 'rRN_inflation', 'polDefAdjusted_inflation', 'rRNwithRiskFree_inflation', 'rRF_inflation'

% Model will be estimated?
conf_estimateModel = false; % true, false

% Include observables even if not estimated
conf_includeObservables = true; % true, false

% modelInterpretation
conf_modelInterpretation = 'Forecasters';

% Turn Off shocks, specify the std parameter (i.e. "sigmaTau")
conf_shocksTurnedOff = ["sigmaTau", "sigmaBeta"];

% List of observables
conf_observablesList = [...%"obs_y", ...
                        ...%"obs_c", ...
                        ...%"obs_g", ...
                        ...%"obs_nr", ...
                        ...%"obs_swap_PreDI_3m", ...
                        ...%"obs_swap_PreDI_6m", ...
                        ...%"obs_unempRate",...
                        ...%"obs_wage", ...
                        ...%"obs_pii", ...
                        ...%"obs_pii_FOCUS_Median", ...
                        ...%"obs_y_FOCUS_Median_BOQ", ...
                        ...%"obs_y_FOCUS_25", ...
                        ...%"obs_y_FOCUS_75", ...
                        ...%"obs_y_FOCUS_0_EOQ", ...
                        ...%"obs_y_FOCUS_20_EOQ", ...
                        ...%"obs_y_FOCUS_25_EOQ", ...
                        ...%"obs_y_FOCUS_50_EOQ", ...
                        ...%"obs_y_FOCUS_75_EOQ", ...
                        ...%"obs_y_FOCUS_80_EOQ", ...
                        ...%"obs_y_FOCUS_100", ...
                        ...%"obs_y_FOCUS_0_EOQ", ...
                        ...%"obs_y_FOCUS_100_EOQ", ...
                        ];
                    
% List of measurement errors
conf_measurementErrorsList = [  
                                ...%"obs_y", ...                              
                                ...%"obs_pii", ...
                                ...%"obs_swap_PreDI_6m"
                                ...%"obs_wage", ...
                                %"obs_y_FOCUS_25", ...
                                %"obs_y_FOCUS_75", ...
                                %"obs_unempRate"
                                %"obs_pii_FOCUS_Median", ...
                                ];                    
% Measurement errors
conf_measurementErrors = true; % true, false

% Shock on effective productivity
conf_shockOnEffectiveProductivity = "None" ; % "None", "OnlyExpectation", "OnlyRealization"

% Shock on productivity observed with delay
conf_psiShockObservedWithDelay = false;

% Build configure struct
riseConfig = struct( ...
        'conf_calibrationName',     conf_calibrationName, ...
        'conf_occBinConstraint',    conf_occBinConstraint, ...
        'conf_regimeSwitching',     conf_regimeSwitching, ...
        'conf_hasGovernment',       conf_hasGovernment, ...
        'conf_govBondsAreRiskFree', conf_govBondsAreRiskFree, ...
        'conf_govExpenses',         conf_govExpenses, ...
        'conf_debtLevel',           conf_debtLevel, ...
        'conf_defPolTarget',        conf_defPolTarget, ...
        'conf_debtLevelGuess',      conf_debtLevelGuess, ...
        'conf_govAccumDebt',        conf_govAccumDebt, ...
        'conf_kVariable',           conf_kVariable, ...
        'conf_prodProcessType',     conf_prodProcessType, ...
        'conf_defRProcessType',     conf_defRProcessType, ...
        'conf_policyRule',          conf_policyRule, ...
        'conf_includeObservables',  conf_includeObservables, ...
        'conf_estimateModel',       conf_estimateModel, ...
        'conf_modelInterpretation', conf_modelInterpretation, ...
        'conf_observablesList',     conf_observablesList, ...
        'conf_measurementErrors',     conf_measurementErrors, ...
        'conf_measurementErrorsList', conf_measurementErrorsList, ...
        'conf_shockOnEffectiveProductivity', conf_shockOnEffectiveProductivity, ...
        'conf_psiShockObservedWithDelay', conf_psiShockObservedWithDelay, ...
        'conf_shocksTurnedOff', conf_shocksTurnedOff, ...
        'conf_simplifyPolicyRule', conf_simplifyPolicyRule,  ...
        'conf_approximateDefaultProb', conf_approximateDefaultProb, ...
        'conf_appendWelfareEquations', conf_appendWelfareEquations, ...
        'conf_stickyPrices', conf_stickyPrices, ...
        'conf_effectiveLowerBound', conf_effectiveLowerBound, ...
        'conf_taxRule', conf_taxRule ...
        );

%% Load data

reloadData = false;

if paperNumber == 2
    startDate   = datetime('01/04/1999', 'InputFormat', 'dd/MM/yyyy'); % one period before first period
    endDate     = datetime('01/10/2019', 'InputFormat', 'dd/MM/yyyy');
else
    startDate   = datetime('01/01/1980', 'InputFormat', 'dd/MM/yyyy');
    endDate     = datetime('01/10/1998', 'InputFormat', 'dd/MM/yyyy'); % one period before first period
end

defLCMeasure = 'EMBI'; % 'CDS', 'EMBI'
dataTransformation = 'diffDemean'; % '', diffDemean, detrendHP

% { {Variable to be demeaned, Variable whose mean is to be used} }
demeanList = {};
% demeanList = {  {'obs_y',       'obs_y'},	...
%                 {'obs_c',       'obs_y'}, ...
%                 {'obs_g',       'obs_y'}, ...
%                 {'obs_inv',     'obs_y'}, ...
%                 {'obs_x',       'obs_y'}, ...
%                 {'obs_imp',     'obs_y'}, ...
%                 };
% { {Variable to be detrended, Variable whose trend is to be used} }
detrendList = {};
% detrendList = { {'obs_y',       'obs_y'},	...
%                 {'obs_c',       'obs_y'}, ...
%                 {'obs_g',       'obs_y'}, ...
%                 {'obs_inv',     'obs_y'}, ...
%                 {'obs_x',       'obs_y'}, ...
%                 {'obs_imp',     'obs_y'}, ...
%                 };

% Pick series to be plotted
if paperNumber == 2
    %obsNames = {'obs_y', 'obs_c', 'obs_g', 'obs_unempRate', 'obs_wage', 'obs_pii', 'obs_defLC', 'obs_swap_PreDI_3m'};
    %obsTexNames = {'$Y^{obs}_t$', '$C^{obs}_t$', '$G^{obs}_t$', 'Unemployment$^{obs}_t$', '$\left(\frac{W_t}{P_t}\right)^{obs}$', '$\pi_t^{obs}$', '$\left(\mathcal{D}_t\right)^{obs}$', '$\left(i_t\right)^{obs}$'};
    %detrendedTexNames = {'$Y^{HP}_t$', '$C^{HP}_t$', '$G^{HP}_t$', 'Unemployment$^{HP}_t$', '$\left(\frac{W_t}{P_t}\right)^{HP}$', '$\pi_t^{HP}$', '$\left(\mathcal{D}_t\right)^{HP}$', '$\left(i_t\right)^{HP}$'};
    obsNames              = {'obs_y', 'obs_c', 'obs_g', 'obs_pii', 'obs_swap_PreDI_3m'};
    obsTexNames           = {'$Y^{obs}_t$', '$C^{obs}_t$', '$G^{obs}_t$', '$\pi_t^{obs}$', '$i_t^{obs}$'};
    beforeDetTexNames     = {'$Y^{obs}_t$', '$C^{obs}_t$', '$G^{obs}_t$', '$\pi_t^{obs}$', '$i_t^{obs}$'};
    detrendedTexNames     = {'$Y^{HP}_t$', '$C^{HP}_t$', '$G^{HP}_t$', '$\pi_t^{deviation}$', '$i_t^{HP}$'};
elseif paperNumber == 3
    obsNames              = {'obs_y', 'obs_pii', 'obs_nr'};
    obsTexNames           = {'$Y^{obs}_t$', '$\pi_t^{obs}$', '$\left(i_t\right)^{obs}$'};
    beforeDetTexNames     = {'$Y^{obs}_t$', '$\pi_t^{obs}$', '$\left(i_t\right)^{obs}$'};
    detrendedTexNames     = {'$Y^{HP}_t$', '$\pi_t^{deviation}$', '$\left(i_t\right)^{HP}$'};
end
plotSelectedVars      = obsNames;
plotSelectedVarsNames = obsTexNames;
plotInterpreter       = 'latex';

if conf_estimateModel && reloadData
    loadData;
else
    load([pathSaved, filesep, 'dataset.mat']);
end

%% Simulate fiscal limits

%addpath('Fiscal Limits');
%distFiscalLimits;

%% Parameterize, solve the model, and assign the data (Policy Rules)

solveOrder  = 1;
solve_derivatives_type = 'symbolic'; % symbolic, automatic, numerical  %%% in case of optimal rule, cannot be symbolic
debugMode   = true; % true, false
steady_state_imposed = false;  % if imposed, then does not check whether it is a steady state (solves around the reference regime ss)
steady_state_unique = false; % the steady state is unique: solves around the "ergodic mean" and does not run the steady state file
solveModel;

%% Goodness of fit of the logistic approximation

f = figure;
previousInterpreter = edu_GraphSetInterpreter('latex');

% Set layout and reduce empty spaces in the graph
tl = tiledlayout(2, 1);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';

custom_LineColor = {'blue', 'red', 'black', 'blue', 'blue', 'blue', 'blue'};

% Annual debt
xRng = 0.8:0.01:1.5;

iReg = 1;
title_defTarget = num2str(mdlVector(1).user_data.conf_defPolTarget);

% Plot interpolated
nexttile(1, [1 1]);
y_p1 = feval(@(b) fFiscalLimit_defProb(4*b,0,0,0,0,0), xRng);
p1 = plot(xRng.*100, y_p1, 'DisplayName', 'Fiscal limits simulation');
p1.Color = 'black';
    
for iMdl = 1:length(mdlVector)    
    
    paramsStructTemp = get(mdlVector(iMdl),'parameters');
    f_defProb = @(b) 1/(1 + exp( paramsStructTemp.probDefFisLim_Param_0(iReg) ...
        + paramsStructTemp.probDefFisLim_Param_B(iReg)*( 4*b - paramsStructTemp.muFisLim_Intercept(iReg)) )) ;
     
    % Plot logistic approximation
    hold on;
    y_p2 = arrayfun(@(x) f_defProb(x), xRng);
    p2 = plot(xRng.*100, y_p2, 'DisplayName', 'Logistic approximation');
    p2.Color = custom_LineColor{iMdl};
    p2.LineStyle = '--';
    p2.LineWidth = 2;
    
    title('Default probability CDF');
    ylabel('cumulative probability');
    xlabel('$\frac{B_t}{4 Y_t}$');
    xlim(100.*[xRng(1) xRng(end)]);
    legend('Location', 'best');
    set(gca, 'FontSize', 12);

    % Plot difference in logistic approximation
    nexttile();
    p3 = plot(xRng.*100, (y_p2 - y_p1) .* 100);
    pub_GraphDrawZeroAxis(p3);
    title('Discrepancy');
    ylabel('probability (\%)');
    xlabel('$\frac{B_t}{4 Y_t}$');
    xlim(100.*[xRng(1) xRng(end)]);
    set(gca, 'FontSize', 12);
end

set(f, 'Position',  [100, 300, 500, 400]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Default Probability' ['Graph - Default Prob.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);
                
%% List parameters

iMdl = 1;
get(mdlVector(iMdl),'parameters')

%% Filter and smooth calibrated model

% Remember that one needs at least as many shocks as observables!

% 1- filtered_variables refer to one-step ahead forecasts: a_{t|t-1}
% 2- updated_variables refer to the updates a_{t|t}
% 3- smoothed_variables refer to the final estimates conditional on all
% available information a_t{t|n}
% 
% As for the expected counterparts, they are averages across all regimes.
% This means that in a regime switching model, you will have for every
% variable, as many series as the number of regimes. The weights used to
% compute those averages are the probabilities.

%close all;

mdlFilter = mdlVector; % mdlVector, estMdl

saveGraphs = false;

firstPlotPeriod = 5;
plotObsAndSmoothed = true;

%%% Create containers of parameters
create_containers = @(n)arrayfun(@(x)containers.Map(), 1:n, 'UniformOutput', false);
mapCalibParams = create_containers(length(mdlVector));
                    
for iMap = 1:length(mdlFilter)
    mdlParams = get(mdlFilter(iMap),'parameters');
    mapCalibParams{iMap} = containers.Map(fieldnames(mdlParams), structfun(@(x) x(1), mdlParams'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlVector = {};
for iMdl = 1:length(mdlFilter)

    % Select variables
    % {'obs_name', 'varName', movAvgPer}   movAvgPer ex: 0 or [2,2]
    allVars             = {
                            'obs_y',                'y',    0; ...
                            'obs_c',                'c',    0; ...
                            'obs_g',                'g',    0; ...
                            ...%'obs_unempRate',...
                            'obs_pii',              'Pii',  0; ...
                            ...'obs_nr', ...
                            'obs_swap_PreDI_3m',    'nrPolicy',     0; ...
                            ...'',                     'rPolicy',      0; ...
                            ...'obs_defLC', 'polDef', 0;...
                            ...'obs_debtOutput_new', 'bY', 0;...
                            ...'obs_tau', 'tau', 0;
                            };
    
    obsVars             = allVars(:,1);
    filterSmoothVars    = allVars(:,2);
    movingAverageVars   = allVars(:,3);
    vlocs               = locate_variables(filterSmoothVars, mdlFilter(iMdl).endogenous.name);
    vtexNames           = mdlFilter(iMdl).endogenous.tex_name(vlocs);
    plotVars            = filterSmoothVars;
    plotVarNames        = vtexNames;

    % Filter and smooth series
    % yRegimes = ones(size(myData,1),1);
    % IMPOSE REGIMES: filter(mdlFilter(iMdl), 'kf_user_algo', {@myKnownRegimesFilter, yRegimes}, 'data', myData)
    if sum(mdlFilter(iMdl).markov_chains.chain_is_endogenous) == 0
        myFiltSmooth = filter(mdlFilter(iMdl), 'data', myData);
    else
        yRegimes = ones(size(myData,1),1);
        %myData2 = myData;
        %myData2.data(:,end) = 0.045 + edu_Detrend(myData2('obs_swap_PreDI_3m').data, 'linear');
        myFiltSmooth = filter(mdlFilter(iMdl), 'kf_user_algo', {@myKnownRegimesFilter, yRegimes}, 'data', myData);
    end

    nCols = 2;
    nLins = ceil(length(plotVars) / nCols);

    f = figure;
    previousInterpreter = edu_GraphSetInterpreter('latex');

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(nLins, nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';


    for iVar = 1:numel(plotVars)

        tl_h = nexttile();
        %tl_h = [tl_h, nexttile()];

        obsData             = myData(obsVars{iVar});
        obsData.varnames    = obsVars(iVar);
        
        if sum(mdlFilter(iMdl).markov_chains.chain_is_endogenous) == 0
            filtData      = myFiltSmooth.filtered_variables.(filterSmoothVars{iVar});
            smoothData    = myFiltSmooth.smoothed_variables.(filterSmoothVars{iVar});
        else
            filtData = myFiltSmooth.filtered_variables.(filterSmoothVars{iVar});
            smoothData = myFiltSmooth.smoothed_variables.(filterSmoothVars{iVar});
            
            filtData.data = filtData.data(1:end-1,1);
            smoothData.data = smoothData.data(:,1);
            
            filtData.varnames = {'filtered'};
            smoothData.varnames = {'smoothed'};
            
%             filtData = NaN(size(dummySmoothData,1),1);
%             smoothData = NaN(size(dummySmoothData,1),1);
%             for iRow = 1:size(dummySmoothData,1)
%                 filtData(iRow)      = dummyFiltData(iRow, yRegimes(iRow));
%                 smoothData(iRow)    = dummySmoothData(iRow, yRegimes(iRow));
%             end
            
        end
        consolData    = [obsData, filtData, smoothData];
        
        seriesRef = consolData;
        dates = (datetime(seriesRef.start, 'InputFormat', 'yyyyQQQ'):calquarters(1):datetime(seriesRef.finish, 'InputFormat', 'yyyyQQQ'))';
        [obsData, filtData, smoothData] = transformVariables(filterSmoothVars{iVar}, mapCalibParams{iMdl}, consolData.data(:,1), consolData.data(:,2), consolData.data(:,3));

        % Moving averages
        % Real interest rate
        if movingAverageVars{iVar} > 0
            filtData = movmean(filtData, [movingAverageVars{iVar}(1) + movingAverageVars{iVar}(2), 0], 'omitnan', 'Endpoints', 'fill');
            smoothData = movmean(smoothData, [movingAverageVars{iVar}(1), movingAverageVars{iVar}(2)], 'omitnan', 'Endpoints', 'fill');
        end
        
        
        nReg = size(filtData, 2);
        plotData = [];
        legStr = {'observed'};
        for iReg = 1:nReg
            plotData = [plotData, real([ ...
                            obsData, ...
                            filtData(:,iReg), ...
                            smoothData(:,iReg) ])
                         ];
            legStr = [legStr {'filtered'} {'smoothed'}];
        end
        graph = plot(dates(firstPlotPeriod:end), plotData(firstPlotPeriod:end,:));
        graph(1).LineWidth = 1;
        graph(2).LineWidth = 1;
        graph(3).LineWidth = 2;
        graph(1).LineStyle = '-';
        graph(2).LineStyle = '--';
        graph(3).LineStyle = ':';
        title(plotVarNames{iVar});
        if iVar == 1
            legend(legStr, 'Location', 'best');
        end
        ylabel('\%');

        edu_GraphDrawZeroAxis(graph);
        
        if plotObsAndSmoothed
            delete(graph(2));
        end
        
        % Row label
        %text(graph(1), mean(xlim(graph(1))), mean(ylim(graph(1))), ...
        %    shocksToPlotNames{iShock}, ...
        %    'Interpreter', 'latex', ...
        %    'FontSize', 20);

        % First row
        %if iVar == 1
        %    graph(1).Legend.Location = 'North';
        %end

        % Hide after the first row
        %if iShock > 1
        %    title(graph(1:end), '');
        %    legend(graph(1), 'off');
        %end

        % Between lines 1 and 3
        %if iShock >= 1 && iShock <= 3
        %    set(graph(3:5), 'Color', [0.9 0.9 0.9]);
        %end

        % Hide after column 2
        %ylabel(graph(3:end), '');

        % Hide after column 2
        %yticks(graph(3:end), []);

    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize legend
    ax = nexttile(1);
    ax.Legend.FontSize = 14;

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 18;
    end

    % Keep axes
    tlVector{iMdl} = findobj(tl.Children, 'Type', 'Axes');
    
    % Fetch config
    modelName   = riseModelNames{iMdl};
    defProcess  = riseConfigVector(iMdl).conf_defRProcessType;
    if riseConfigVector(iMdl).conf_govBondsAreRiskFree
        riskyPolicy = 'Risk-Free';
    else
        riskyPolicy = 'Risky';
    end

    set(f, 'Position',  [100, 100, 800, 800]); % resize figure
    if saveGraphs
        exportgraphics(f, ...
            strjoin({pathImages 'Filtered' ['Graph - Filtered and Smoothed Shocks'  ' - DebtLevel ' debtLevel '.png']}, filesep));
    end
    edu_GraphSetInterpreter(previousInterpreter);
    
end

% Link axes
if length(mdlFilter) > 1
    for iVar = 1:numel(plotVars)
        arrAxes = [];
        for iMdl = 1:length(mdlFilter)
            arrAxes = [arrAxes, tlVector{iMdl}(iVar)];
        end
        linkaxes(arrAxes);
    end
end

%% Plot rPolicy

% 1- filtered_variables refer to one-step ahead forecasts: a_{t|t-1}
% 2- updated_variables refer to the updates a_{t|t}
% 3- smoothed_variables refer to the final estimates conditional on all
% available information a_t{t|n}
% 
% As for the expected counterparts, they are averages across all regimes.
% This means that in a regime switching model, you will have for every
% variable, as many series as the number of regimes. The weights used to
% compute those averages are the probabilities.

close all;

firstPlotPeriod = 5;
plotObsAndSmoothed = true;

%%% Create containers of parameters
create_containers = @(n)arrayfun(@(x)containers.Map(), 1:n, 'UniformOutput', false);
mapCalibParams = create_containers(length(mdlVector));
                    
for iMap = 1:length(mdlVector)
    mdlParams = get(mdlVector(iMap),'parameters');
    mapCalibParams{iMap} = containers.Map(fieldnames(mdlParams), structfun(@(x) x(1), mdlParams'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select variables
if paperNumber == 2
    obsVars             = {'obs_y', 'obs_swap_PreDI_3m', 'obs_r', 'obs_defLC'};
    filterSmoothVars    = {'y', 'nrPolicy', 'rPolicy', 'polDef'};
    movingAverageVars    = {[1,1], [1,1], [1,1], [1,1]};
elseif paperNumber == 3
    obsVars             = {'obs_y', 'obs_nr', 'obs_r', 'obs_defLC'};
    filterSmoothVars    = {'y', 'nrPolicy', 'rPolicy', 'polDef'};
    movingAverageVars    = {[1,1], [1,1], [1,1], [1,1]};
end

tlVector = {};
for iVar = 1:numel(obsVars)
    
    f = figure;
    previousInterpreter = edu_GraphSetInterpreter('latex');
    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(1, 1);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';
    
    legStr = {};
    for iMdl = 1:length(mdlVector)

        vlocs               = locate_variables(filterSmoothVars, mdlVector(iMdl).endogenous.name);
        vtexNames           = mdlVector(iMdl).endogenous.tex_name(vlocs);
        plotVars            = filterSmoothVars;
        plotVarNames        = vtexNames;

        % Filter and smooth series
        myFiltSmooth = filter(mdlVector(iMdl),'data',data);
    
        %tl_h = nexttile();
        %tl_h = [tl_h, nexttile()];

        obsData       = data.(obsVars{iVar});
        filtData      = myFiltSmooth.filtered_variables.(filterSmoothVars{iVar});
        smoothData    = myFiltSmooth.smoothed_variables.(filterSmoothVars{iVar});
        consolData    = [obsData, filtData, smoothData]; 

        seriesRef = consolData;
        dates = (datetime(seriesRef.start, 'InputFormat', 'yyyyQQQ'):calquarters(1):datetime(seriesRef.finish, 'InputFormat', 'yyyyQQQ'))';
        [obsData, filtData, smoothData] = transformVariables(filterSmoothVars{iVar}, mapCalibParams{iMdl}, consolData.data(:,1), consolData.data(:,2), consolData.data(:,3));

        % Moving averages
        % Real interest rate
        if movingAverageVars{iVar} > 0
            filtData = movmean(filtData, [movingAverageVars{iVar}(1) + movingAverageVars{iVar}(2), 0], 'omitnan', 'Endpoints', 'fill');
            smoothData = movmean(smoothData, [movingAverageVars{iVar}(1), movingAverageVars{iVar}(2)], 'omitnan', 'Endpoints', 'fill');
        end
        
        
        nReg = size(filtData, 2);
        plotData = [];
        if iMdl == 1
            legStr = [legStr, {'observed'}];
            
            for iReg = 1:nReg
                plotData = [plotData, real([ ...
                                obsData, ...
                                filtData(:,iReg), ...
                                smoothData(:,iReg) ])
                             ];
                legStr = [
                    legStr , ...
                    ...%{[riseSetUpNames{iMdl}, ' - filtered']}, ...
                    {[riseSetUpNames{iMdl}, ' - smoothed']}
                    ];
            end
            hold on;
            graph = plot(dates(firstPlotPeriod:end), plotData(firstPlotPeriod:end,:));
            hold off;
            graph(1).LineWidth = 1;
            graph(2).LineWidth = 1;
            graph(3).LineWidth = 2;
            graph(1).LineStyle = '-';
            graph(2).LineStyle = '--';
            graph(3).LineStyle = ':';
            
            if plotObsAndSmoothed
                delete(graph(2));
            end
        else
            for iReg = 1:nReg
                plotData = [plotData, real([ ...
                                filtData(:,iReg), ...
                                smoothData(:,iReg) ])
                             ];
                legStr = [
                    legStr, ...
                    ...%{[riseSetUpNames{iMdl}, ' - filtered']}, ...
                    {[riseSetUpNames{iMdl}, ' - smoothed']}];
            end
            hold on;
            graph = plot(dates(firstPlotPeriod:end), plotData(firstPlotPeriod:end,:));
            hold off;
            
            graph(1).LineWidth = 1;
            graph(2).LineWidth = 2;
            graph(1).LineStyle = '--';
            graph(2).LineStyle = ':';
            
            if plotObsAndSmoothed
                delete(graph(1));
            end
            
        end

    end

    edu_GraphDrawZeroAxis(gca);
    
    title(plotVarNames{iVar});
    legend(legStr, 'Location', 'best');
    ylabel('\%');
    
%     % Custom font
%     set(f.Children.Children, 'FontSize', 12);
%     set(f.Children.Children, 'FontWeight', 'bold');

    % Resize legend
    %ax.Legend.FontSize = 14;
    
    % Resize title
    for iGraph = 1:length(tl.Children)
        graph.Parent.Title.FontSize = 18;
        set(graph.Parent.Children, 'FontSize', 18);
        set(graph.Parent.Children, 'FontWeight', 'bold');
    end

    % Fetch config
    modelName   = riseModelNames{1};
    defProcess  = riseConfigVector(1).conf_defRProcessType;

    set(f, 'Position',  [100, 100, 800, 800]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'IRFs' modelName 'DefR' defProcess ['Graph - All Rules - ' filterSmoothVars{iVar} '.png']}, filesep));
    edu_GraphSetInterpreter(previousInterpreter);

end

%% Estimate model

estimateModel;

%% Estimation summary

estimParamNames = fieldnames(estimPostSimData);
for ii = 1:length(estimParamNames)
    pStruct = estimPostSimData.(estimParamNames{ii});
    disp( [ estimParamNames{ii} ' 5%  '  num2str(pStruct.x_kdens(find(cumsum(pStruct.f_kdens)/sum(pStruct.f_kdens) >= 0.05, 1))) ...
            ' 95%  ' num2str(pStruct.x_kdens(find(cumsum(pStruct.f_kdens)/sum(pStruct.f_kdens) >= 0.95, 1)))
        ]);
end

%% Print solution

for iMdl = 1:nMdl
    mdlVector(iMdl).print_solution();
end    

%% Solution residuals per equation

residMatrix     = [];
residVarNames   = [];
for iMdl = 1:nMdl
    nRegimes    = mdlVector(iMdl).markov_chains.regimes_number;
    nEquations  = mdlVector(iMdl).endogenous.number;

    residMatrix     = [residMatrix, resid(mdlVector(iMdl))];
    residVarNames   = [residVarNames; strcat("Mdl ", num2str(iMdl), "; Regime ", string(1:nRegimes)')];
end

tRes = array2table(residMatrix);
tRes.Properties.VariableNames = residVarNames;
tRes(:,'Equation') = mdlVector(1).equations.dynamic;
tRes = tRes(:,["Equation"; residVarNames]);
tRes.Properties.RowNames = strcat("Equation ", string(1:nEquations)');
disp(tRes);

%% Steady state

mdlVector(1).markov_chains.regimes

ssMatrix     = [];
ssVarNames   = [];
for iMdl = 1:nMdl
    nRegimes     = mdlVector(iMdl).markov_chains.regimes_number;
    ssMatrix     = [ssMatrix, [mdlVector(iMdl).solution.ss{:}]];
    ssVarNames   = [ssVarNames; strcat("Mdl ", num2str(iMdl), "; Regime ", string(1:nRegimes)')];
end
tSS = array2table(ssMatrix);
tSS.Properties.VariableNames    = ssVarNames;
tSS.Properties.RowNames         = mdlVector(1).endogenous.name;
disp(tSS);

%% Evaluate solution existence and stability

% Is there a solution?
disp('Are there solutions?');
nSols = NaN(nMdl,1);
for ii = 1:nMdl
    nSols(ii) = mdlVector(ii).nsols;
end
tSolutions = table(nSols, ...
    'VariableNames', {'nSolutions'}, 'RowNames', strcat("Model ", string(1:length(riseConfigVector))'));
disp(tSolutions)

% Is stationary?
disp('Is the linear markov switching system stationary?');
stab = NaN(nMdl,1);
for ii = 1:nMdl
    stab(ii) = is_stationary_system(mdlVector(ii));
end
tStationarity = table(stab, ...
    'VariableNames', {'isStable'}, 'RowNames', strcat("Model ", string(1:length(riseConfigVector))'));
disp(tStationarity)

% Is stable?
% Checks the stability of a linear markov switching system
disp('Is the linear markov switching system stable?');
cfm = NaN(nMdl,1);
gmh = NaN(nMdl,1);
for ii = 1:nMdl
    stability_algorithm = 'cfm'; % 'cfm' or 'gmh'
    cfm(ii) = is_stable_system(mdlVector(ii), 'stability_algorithm', stability_algorithm);
    stability_algorithm = 'gmh'; % 'cfm' or 'gmh'
    gmh(ii) = is_stable_system(mdlVector(ii), 'stability_algorithm', stability_algorithm);
end
tStability = table(cfm, gmh, ...
    'VariableNames', {'cfm', 'gmh'}, 'RowNames', strcat("Model ", string(1:length(riseConfigVector))'));
disp(tStability)

%% Test stability

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlStabVector = mdlVector;
for iMdl = 1:length(mdlStabVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlStabVector(iMdl)      = set(mdlStabVector(iMdl),'parameters', paramsStructTemp);
    mdlStabVector(iMdl)      = solve(mdlStabVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

% Load parameters / Simulate fiscal limit
loadParamsAndSimulFiscalLimits;

% Override calibration parameters
paramNamesOverride      = {};
paramValuesOverride     = {};

for iMdl = 1:length(mdlStabVector)
    
    config  = riseConfigVector(iMdl);
    paramNames   = {'phi', 'gammaTau'};
    %paramNames   = {'phi', 'gammaTau'};
    paramPlotNames = {'$\phi^{\pi}$', '$\gamma_\tau$'};
    paramValues = { 0:0.1:3.0, 0:0.01:0.30};
    stability_algorithm = 'cfm'; % 'cfm' or 'gmh'
    stateValues     = {-1,0,1,2,3};
    stateNames      = {'No solution', 'Solution is unstable', 'Solution is stable', 'Multiple solutions', 'No solution found'};
    
    stabResults = NaN(length(paramValues{1}), length(paramValues{2}));
    for iParam = 1:length(paramValues{1})
        for jParam = 1:length(paramValues{2})
            disp(['Checking stability with parameter value = ', ...
                num2str(paramValues{1}(iParam)) ' and ' num2str(paramValues{2}(jParam))]);
            
            paramsStructTemp = paramsStruct;
            paramsStructTemp.(paramNames{1}) = paramValues{1}(iParam);
            paramsStructTemp.(paramNames{2}) = paramValues{2}(jParam);

            % Override parameters
            for iParamOver = 1:length(paramNamesOverride)
                paramsStructTemp.(paramNamesOverride{iParamOver}) = paramValuesOverride{iParamOver};
            end
            
            mdl     = mdlStabVector(iMdl);
            mdl     = set(mdl,'parameters', paramsStructTemp);
            eigValsStr = evalc('mdl = solve(mdl);');
            disp(eigValsStr);
            
            if mdl.nsols > 0
                if is_stable_system(mdl, 'stability_algorithm', stability_algorithm)
                     % Stable solution
                    stabResults(iParam, jParam) = stateValues{3};
                else
                    % Non-stable solution
                    stabResults(iParam, jParam) = stateValues{2};
                end
            else
                eigvals = edu_ExtractNumberFromString(eigValsStr);
                if isempty(eigvals)
                    % Regime-switching with no solution found
                    stabResults(iParam, jParam) = stateValues{5};
                elseif ~isempty(eigvals) && eigvals(1) > eigvals(3)
                    % Multiple solutions
                    stabResults(iParam, jParam) = stateValues{4};
                else
                    % No solution
                    stabResults(iParam, jParam) = stateValues{1};
                end
            end
        end
    end

    % Filter used states
    usedStateValues = unique(stabResults);
    usedStateNames  = stateNames(find(ismember(cell2mat(stateValues), usedStateValues)));

    % Plot parameter stability
    xN = length(paramValues{1});
    yN = length(paramValues{2});
    iItem = 0;
    dMatrix = NaN(xN*yN, 3);
    for ii = 1:length(paramValues{1})
        for jj = 1:length(paramValues{2})
            iItem = iItem + 1;
            dMatrix(iItem, 1) = paramValues{1}(ii);
            dMatrix(iItem, 2) = paramValues{2}(jj);
            dMatrix(iItem, 3) = stabResults(ii, jj);
        end
    end

    f = figure;
    previousInterpreter = edu_GraphSetInterpreter('latex');
    h = gscatter(dMatrix(:,2), dMatrix(:,1), dMatrix(:,3));
    for iGroup = 1:length(unique(dMatrix(:,3)))
        h(iGroup).Marker = 's';
        h(iGroup).MarkerSize = 5;
        h(iGroup).MarkerEdgeColor = [0 0 0];
        switch h(iGroup).DisplayName
            case '-1'
                h(iGroup).Marker = 'o';
                h(iGroup).MarkerFaceColor = [0.7, 0.1, 0.0];
                h(iGroup).Color = [0.7, 0.1, 0.0];
            case '0'
                h(iGroup).Marker = '.';
                h(iGroup).MarkerFaceColor = [1.0, 0.9, 0.1];
                h(iGroup).Color = [1.0, 0.9, 0.1];
            case '1'
                h(iGroup).Marker = '+';
                h(iGroup).MarkerFaceColor = [0.0, 0.3, 0.4];
                h(iGroup).Color = [0.0, 0.3, 0.4];
            case '2'
                h(iGroup).Marker = '*';
                h(iGroup).MarkerFaceColor = [0.0, 0.0, 0.0];
                h(iGroup).Color = [0.0, 0.0, 0.0];
            case '3'
                h(iGroup).Marker = 'o';
                h(iGroup).MarkerFaceColor = [0.7, 0.1, 0.0];
                h(iGroup).Color = [0.7, 0.1, 0.0];
        end
    end
    xlabel(paramPlotNames{2});
    ylabel(paramPlotNames{1});
    xticks(paramValues{2}(1):2*(paramValues{2}(2) - paramValues{2}(1)):paramValues{2}(end));
    yticks(paramValues{1}(1):2*(paramValues{1}(2) - paramValues{1}(1)):paramValues{1}(end));
    xlim([
        paramValues{2}(1) - (paramValues{2}(2) - paramValues{2}(1))/2, ...
        paramValues{2}(end) + (paramValues{2}(2) - paramValues{2}(1))/2
        ]);
    ylim([
        paramValues{1}(1) - (paramValues{1}(2) - paramValues{1}(1))/2, ...
        paramValues{1}(end) + (paramValues{1}(2) - paramValues{1}(1))/2
        ]);
    legend(usedStateNames, 'Location', 'northoutside', 'Orientation', 'horizontal');
    set(gca, 'TickLength',[0 0])
    set(gca, 'FontSize', 14);
    edu_GraphSetInterpreter(previousInterpreter);
    
    simGraphName = ['ParamStability_' paramNames{1} '_and_' paramNames{2} '_' mdlStabVector(iMdl).user_data.conf_policyRule ' - DefProb ' num2str(mdlStabVector(iMdl).user_data.conf_defPolTarget) '.png'];
  
    set(gcf, 'Position',  [100, 100, 800, 600]); % resize figure
    exportgraphics(f, ...
    strjoin({pathImages 'ParamStability' title_stickyPrices simGraphName}, filesep));

end

%% Test stability and sign identification

for ii = 1:1
    
    mdl     = mdlVector(ii);
    config  = riseConfigVector(ii);
    paramNames   = {'sigma', 'eta'};
    paramPlotNames = {'$\sigma$', '$\eta$'};
    paramValues = {0.8:0.1:1.5, 1.0:2.0:11.0};
    stability_algorithm = 'cfm'; % 'cfm' or 'gmh'
    stateValues = {-1,0,1};
    stateNames  = {'No equilibrium', 'Equilibrium is unstable', 'Equilibrium is stable'};

    stabResults         = NaN(length(paramValues{1}), length(paramValues{2}));
    rRNResults          = NaN(length(paramValues{1}), length(paramValues{2}), length(shock_list));
    rGoodResults        = NaN(length(paramValues{1}), length(paramValues{2}), length(shock_list));
    rBadResults         = NaN(length(paramValues{1}), length(paramValues{2}), length(shock_list));
    rGoodOnlyResults    = NaN(length(paramValues{1}), length(paramValues{2}), length(shock_list));
    for iParam = 1:length(paramValues{1})
        for jParam = 1:length(paramValues{2})
            disp(['Checking stability with parameter value = ', ...
                num2str(paramValues{1}(iParam)) ' and ' num2str(paramValues{2}(jParam))]);
            paramsStruct.(paramNames{1}) = paramValues{1}(iParam);
            paramsStruct.(paramNames{2}) = paramValues{2}(jParam);

            mdl = set(mdl,'parameters', paramsStruct);
            mdl = solve(mdl);

            if mdl.nsols > 0
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                shockSize = 1;
                M = set(mdl,  'irf_shock_sign', shockSize);

                nIRFPeriods                             = 100;
                irf_type                                = 'irf'; % 'irf', 'girf'
                simul_honor_constraints                 = true;
                simul_honor_constraints_through_switch  = true;
                irf_regime_specific                     = false;
                
                computeIRFs;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for iShock = 1:length(shock_list)
                    rRNResults(iParam, jParam, iShock)          = sign(myirfs.(shock_list{iShock}).rRN(1,1).values);
                    rGoodResults(iParam, jParam, iShock)        = sign(myirfs.(shock_list{iShock}).rGood(1,1).values);
                    rBadResults(iParam, jParam, iShock)         = sign(myirfs.(shock_list{iShock}).rBad(1,1).values);
                    rGoodOnlyResults(iParam, jParam, iShock)    = sign(myirfs.(shock_list{iShock}).rGoodOnly(1,1).values);
                end
                
                if is_stable_system(mdl, 'stability_algorithm', stability_algorithm)
                    stabResults(iParam, jParam) = stateValues{3};
                else
                    stabResults(iParam, jParam) = stateValues{2};
                end
            else
                stabResults(iParam, jParam) = stateValues{1};
            end
        end
    end

    % Filter used states
    usedStateValues = unique(stabResults);
    usedStateNames  = stateNames(find(ismember(cell2mat(stateValues), usedStateValues)));

    %%%%%%%% Plot parameter stability
    xN = length(paramValues{1});
    yN = length(paramValues{2});
    iItem = 0;
    dMatrix = NaN(xN*yN, 3);
    for ii = 1:length(paramValues{1})
        for jj = 1:length(paramValues{2})
            iItem = iItem + 1;
            dMatrix(iItem, 1) = paramValues{1}(ii);
            dMatrix(iItem, 2) = paramValues{2}(jj);
            dMatrix(iItem, 3) = stabResults(ii, jj);
        end
    end

    f = figure;
    previousInterpreter = edu_GraphSetInterpreter('latex');
    h = gscatter(dMatrix(:,2), dMatrix(:,1), dMatrix(:,3));
    for iGroup = 1:length(unique(dMatrix(:,3)))
        h(iGroup).Marker = 's';
        h(iGroup).MarkerSize = 10;
        h(iGroup).MarkerEdgeColor = [0 0 0];
        switch h(iGroup).DisplayName
            case '-1'
                h(iGroup).MarkerFaceColor = [0 0 0];
            case '0'
                h(iGroup).MarkerFaceColor = [0.5 0.5 0.5];
            case '1'
                h(iGroup).MarkerFaceColor = [1 1 1];
        end
    end
    xlabel(paramPlotNames{2});
    ylabel(paramPlotNames{1});
    xticks(paramValues{2}(1):0.02:paramValues{2}(end));
    yticks(paramValues{1}(1):0.1:paramValues{1}(end));
    xlim([
        paramValues{2}(1) - (paramValues{2}(2) - paramValues{2}(1))/2, ...
        paramValues{2}(end) + (paramValues{2}(2) - paramValues{2}(1))/2
        ]);
    ylim([
        paramValues{1}(1) - (paramValues{1}(2) - paramValues{1}(1))/2, ...
        paramValues{1}(end) + (paramValues{1}(2) - paramValues{1}(1))/2
        ]);
    legend(usedStateNames, 'Location', 'northoutside');
    set(gca, 'FontSize', 14);
    edu_GraphSetInterpreter(previousInterpreter);

    if config.conf_hasGovernment
            simGraphName = ['withGov_paramStability_' paramNames{1} '_and_' paramNames{2} '_' config.conf_policyRule '.png'];
    else
            simGraphName = ['withoutGov_paramStability_' paramNames{1} '_and_' paramNames{2} '_' config.conf_policyRule '.png'];
    end
    set(gcf, 'Position',  [100, 100, 1000, 1000]); % resize figure
    exportgraphics(f, [pathImages, filesep , 'ParamStability', filesep, simGraphName]);
    %%%%%%%%
    
end

% Show sign basis
[rGoodResults(:,:,1) rBadResults(:,:,1) rGoodOnlyResults(:,:,1)]
[rGoodResults(:,:,2) rBadResults(:,:,2) rGoodOnlyResults(:,:,2)]
[rGoodResults(:,:,3) rBadResults(:,:,3) rGoodOnlyResults(:,:,3)]

%% Simulate Endogenous Regime-Switching

mdl = mdlVector(1);

[mySims,simStates,retcode] = simulate(mdl, ...
        'simul_regime',             simul_regime, ...
        'simul_order',              simul_order, ...
        'simul_pruned',             simul_pruned, ...
        'simul_honor_constraints',  simul_honor_constraints, ...
        'simul_periods',            simul_periods, ...
        'simul_burn',               simul_burn, ...
        'simul_frwrd_back_shoot',   simul_frwrd_back_shoot, ...
        'simul_shock_uncertainty',  simul_shock_uncertainty);

myirfs = irf(mdl, ...
                'irf_periods', nIRFPeriods, ...
                'irf_type', irf_type, ...
                'irf_regime_specific', irf_regime_specific, ...
                'irf_draws', irf_draws, ...
                'irf_girf_regime_uncertainty', irf_girf_regime_uncertainty, ...
                'simul_honor_constraints', simul_honor_constraints, ...
                'simul_honor_constraints_through_switch', simul_honor_constraints_through_switch);

%% Simulate Welfare

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlSimVector = mdlVector;
for iMdl = 1:length(mdlSimVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlSimVector(iMdl)      = set(mdlSimVector(iMdl),'parameters', paramsStructTemp);
    mdlSimVector(iMdl)      = solve(mdlSimVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 4;
simul_regime                            = 1:nRegimes;
simul_order                             = 2;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 30000;
simul_burn                              = 3000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% Load parameters / Simulate fiscal limit
%loadParamsAndSimulFiscalLimits;

% Override calibration parameters
paramNamesOverride      = {};
paramValuesOverride     = {};

% If param is switching, define as i.e. {{'phi_fisLim_1', 'phi_fisLim_2'}} to
% match params in every iteration
paramNames      = {'phi', 'gammaTau'};
paramRef        = {paramsStruct.phi, paramsStruct.gammaTau};
paramPlotNames  = {'\phi', '\gamma_\tau'};
paramTexPlotNames  = {'$\phi^{\pi}$', '$\gamma_\tau$'};
param1Ref = simStruct.params.(paramNames{1});
param2Ref = simStruct.params.(paramNames{2});
paramValues     = { 1.0:0.25:3.0, 0.1:0.025:0.2};
stability_algorithm = 'cfm'; % 'cfm' or 'gmh'

simulateWelfare;

%% Simulate Welfare (Ricardian vs. Non-Ricardian)

nSims                   = 1;
simul_regime            = 1:nRegimes;
simul_order             = 2;
simul_pruned            = true;
simul_honor_constraints = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods           = 1000;
simul_burn              = 200;
simul_frwrd_back_shoot  = false; % true or false
simul_shock_uncertainty = true; % true or false

% Load parameters / Simulate fiscal limit
loadParamsAndSimulFiscalLimits;

% Override calibration parameters
paramNamesOverride      = {};
paramValuesOverride     = {};

paramNames   = {'phi_fisLim_1', 'gammaTau_fisLim_1', 'fracNR'};
%paramNames   = {'phi', 'gammaTau'};
paramRef        = {paramsStruct.phi_fisLim_1, paramsStruct.gammaTau_fisLim_1, paramsStruct.fracNR};
paramPlotNames  = {'\phi', '\gamma_\tau', '\gamma^{NR}'};
paramValues     = { 1.0:0.2:2.0, 0:0.02:0.20, 0:0.1:1.0};
stability_algorithm = 'cfm'; % 'cfm' or 'gmh'

simulateWelfare;

%% Simulate Interest Rates and Inflation

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlSimVector = mdlVector;
for iMdl = 1:length(mdlSimVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlSimVector(iMdl)      = set(mdlSimVector(iMdl),'parameters', paramsStructTemp);
    mdlSimVector(iMdl)      = solve(mdlSimVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 4;
simul_regime                            = 1:nRegimes; % 1:nRegimes
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 30000;
simul_burn                              = 3000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% Load parameters / Simulate fiscal limit
%loadParamsAndSimulFiscalLimits;

% Override calibration parameters
paramNamesOverride      = {}; % ex: tfpLoss, deltaBar
paramValuesOverride     = {}; % ex: 0

% If param is switching, define as i.e. {{'phi_fisLim_1', 'phi_fisLim_2'}} to
% match params in every iteration
regimesToCalculateMoments = 1; % ex: 1:nRegimes vector with regimes to include in the calculation of moments
paramNames      = {'phi', 'phiNr', 'phi_Y'};
paramRef        = {paramsStruct.phi, paramsStruct.phiNr, paramsStruct.phi_Y};
paramPlotNames  = {'\phi', '\phi^i', '\phi^Y'};
paramTexPlotNames  = {'$\phi^{\pi}$', '$\phi^i$', '$\phi^Y$'};
param1Ref = simStruct.params.(paramNames{1});
param2Ref = simStruct.params.(paramNames{2});
param3Ref = simStruct.params.(paramNames{3});
paramValues     = { 1.0:0.25:3.0, 0.0:0.15:0.9, 0.0};
stability_algorithm = 'cfm'; % 'cfm' or 'gmh'

simulateNrPii;

%% How to get the same simulation twice
% Inspiration: https://github.com/jmaih/RISE_toolbox/issues/52

nosimul=3000;
mysims=simulate(mdlVector,'simul_honor_constraints',true,'simul_periods',nosimul,'simul_burn',0); %with 'simul_burn',0, mysims will include the steady state (i.e. all shocks are zero) in the first period

db = initial_conditions(mdlVector, 1, 'steady', 1);

% set initial conditions for shocks
db.epsA(1, 1, 1:(nosimul+1))  = mysims.epsA.data;
db.epsG(1, 1, 1:(nosimul+1))  = mysims.epsG.data;
db.epsM(1, 1, 1:(nosimul+1))  = mysims.epsM.data;

% specify active shocks
shocks = {'epsA','epsG','epsM'};

% set initial conditions for regime
% create the field db.regime; by copying db.e_r to it it will have the right
% character (ts); if the first step is omitted, db.regime will have the
% character double and the next simulation does not run
db.regime = db.epsA;
db.regime(1, 1, 1:(nosimul+1)) = mysims.regime.data;

% Note that due to the usage of "'simul_historical_data',db" and 
% the definition of db as "db = initial_conditions(m, 1, 'steady', 1)" 
% the value of all endogenous variables in the first period will be
% its steady state value
 mysims2=simulate(mdlVector,'simul_honor_constraints',true,'simul_historical_data',db,'simul_periods',nosimul,'forecast_cond_exo_vars',shocks,'simul_shock_uncertainty',false);

 f = figure;
 plot([mysims.c, mysims2.c]);
 
 %% Calculate the Risky Steady State (Under Construction)
% Inspiration: https://github.com/jmaih/RISE_toolbox/issues/52

nSims                                   = 1;
simul_regime                            = 1;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 10000;
simul_burn                              = 1000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false


db = initial_conditions(mdlVector, 1, 'steady', 1);

% set initial conditions for shocks
db.epsA(1, 1, 1:(simul_burn+1))  = 0;
db.epsG(1, 1, 1:(simul_burn+1))  = 0;
db.epsM(1, 1, 1:(simul_burn+1))  = 0;

% specify active shocks
shocks = {'epsA','epsG','epsM'};

% set initial conditions for regime
% create the field db.regime; by copying db.e_r to it it will have the right
% character (ts); if the first step is omitted, db.regime will have the
% character double and the next simulation does not run
%db.regime = db.epsA;
%db.regime(1, 1, 1:(simul_burn+1)) = regIdx;

% Note that due to the usage of "'simul_historical_data',db" and 
% the definition of db as "db = initial_conditions(m, 1, 'steady', 1)" 
% the value of all endogenous variables in the first period will be
% its steady state value
simul_sig = 0;
mysims_0 = simulate(mdlVector, ...
                        'simul_order', 1, ...
                        'simul_honor_constraints', false, ...
                        'simul_historical_data', db, ...
                        'simul_periods', simul_burn, ...
                        'forecast_cond_exo_vars', shocks, ...
                        'simul_shock_uncertainty', true, ...
                        'simul_sig', simul_sig ...
                        );

mysims_1 = simulate(mdlVector, ...
                        'simul_order', 1, ...
                        'simul_honor_constraints', false, ...
                        'simul_historical_data', db, ...
                        'simul_periods', simul_burn, ...
                        'forecast_cond_exo_vars', shocks, ...
                        'simul_shock_uncertainty', true, ...
                        'simul_anticipate_zero', true, ...
                        'simul_honor_constraints_through_switch', true ...
                        );

mysims_2 = simulate(mdlVector, ...
                        'simul_order', 2, ...
                        'simul_honor_constraints', false, ...
                        'simul_historical_data', db, ...
                        'simul_periods', simul_burn, ...
                        'forecast_cond_exo_vars', shocks, ...
                        'simul_shock_uncertainty', true, ...
                        'simul_anticipate_zero', true, ...
                        'simul_honor_constraints_through_switch', true ...
                        );

mysims_3 = simulate(mdlVector, ...
                        'simul_order', 3, ...
                        'simul_honor_constraints', false, ...
                        'simul_historical_data', db, ...
                        'simul_periods', simul_burn, ...
                        'forecast_cond_exo_vars', shocks, ...
                        'simul_shock_uncertainty', true, ...
                        'simul_anticipate_zero', true, ...
                        'simul_honor_constraints_through_switch', true ...
                        );                    
                    
%                         'simul_regime',             simul_regime, ...
%                         'simul_order',              simul_order, ...
%                         'simul_pruned',             simul_pruned, ...
%                         'simul_honor_constraints',  simul_honor_constraints, ...
%                         'simul_honor_constraints_through_switch', simul_honor_constraints_through_switch, ...
%                         'simul_anticipate_zero',    simul_anticipate_zero, ...
%                         'simul_periods',            simul_periods, ...
%                         'simul_burn',               simul_burn, ...
%                         'simul_frwrd_back_shoot',   simul_frwrd_back_shoot, ...
%                         'simul_shock_uncertainty',  simul_shock_uncertainty ...                    
                    
f = figure;
tl = tiledlayout(2,2);
nexttile();
plot(1:simul_burn+1, [mysims_0.c, mysims_1.c, mysims_2.c, mysims_3.c]);
title('C');
nexttile();
plot(1:simul_burn+1, [mysims_0.rPolicy, mysims_1.rPolicy, mysims_2.rPolicy]);
title('rPolicy');
nexttile();
plot(1:simul_burn+1, [mysims_0.nrPolicy, mysims_1.nrPolicy, mysims_2.nrPolicy]);
title('nrPolicy');
nexttile();
plot(1:simul_burn+1, [mysims_0.Pii, mysims_1.Pii, mysims_2.Pii]);
title('Pii');

%% Simulate Correlation: inflation and default prob.

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlCorrVector = mdlVector;
for iMdl = 1:length(mdlCorrVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlCorrVector(iMdl)      = set(mdlCorrVector(iMdl),'parameters', paramsStructTemp);
    mdlCorrVector(iMdl)      = solve(mdlCorrVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 4;
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = false;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 30000;
simul_burn                              = 3000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% Load parameters / Simulate fiscal limit
%loadParamsAndSimulFiscalLimits;

% Override calibration parameters
paramNamesOverride      = {};
paramValuesOverride     = {};

% If param is switching, define as i.e. {{'phi_fisLim_1', 'phi_fisLim_2'}} to
% match params in every iteration
paramNames      = {'phi', 'gammaTau'};
paramRef        = {paramsStruct.phi, paramsStruct.gammaTau};
paramPlotNames  = {'\phi^\pi', '\gamma_\tau'};
paramValues     = { 0.0:0.25:3.0, 0:0.025:0.30};
stability_algorithm = 'cfm'; % 'cfm' or 'gmh'

correlVarNames  = {'Pii', 'fisLim_1_2'};
correlTexVars = {'$\Pi_t$', '$\mathcal{D}_{t+1}$'};

tic;
simulateCorrelation;
toc;

%simulate_rGap;

%% Simulate Probability of Regime-Switching

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlSimVector = mdlVector;
for iMdl = 1:length(mdlSimVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlSimVector(iMdl)      = set(mdlSimVector(iMdl),'parameters', paramsStructTemp);
    mdlSimVector(iMdl)      = solve(mdlSimVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 4;
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 50000;
simul_burn                              = 5000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% Load parameters / Simulate fiscal limit
%loadParamsAndSimulFiscalLimits;

% Override calibration parameters
paramNamesOverride      = {};
paramValuesOverride     = {};

% If param is switching, define as i.e. {{'phi_fisLim_1', 'phi_fisLim_2'}} to
% match params in every iteration
paramNames      = {'phi', 'gammaTau'};
paramRef        = {paramsStruct.phi, paramsStruct.gammaTau};
paramPlotNames  = {'$\phi$', '$\gamma_\tau$'};
paramValues     = { 1.0:0.5:5.0, 0.10:0.05:0.25};
stability_algorithm = 'cfm'; % 'cfm' or 'gmh'

simulateRSProb;

%simulate_rGap;

%% Simulate Probability that Constraints Bind

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 1;
simul_regime                            = 1:4;
simul_order                             = 1;
simul_pruned                            = false;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 500;
simul_burn                              = 100;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% Turn off specific shocks
paramsStructTemp = paramsStruct;
paramsStructTemp.('sigmaM')         = 0;
paramsStructTemp.('sigmaTau')       = 0;
paramsStructTemp.('sigmaPolDef')    = 0;
paramsStructTemp.('sigmaBeta')      = 0;
%paramsStructTemp.('sigmaA')         = 0;
%paramsStructTemp.('sigmaG')         = 0;
mdl = mdlVector;
nMdl = length(mdl);
mdl = set(mdl,'parameters', paramsStructTemp);
mdl = solve(mdl);

% Simulate the model
simulateModel;

% Save sim results
%save([pathSaved filesep 'sim_ProbConstraintsBind.mat']);
%load([pathSaved filesep 'sim_ProbConstraintsBind.mat']);

% Verify regime binding
nStates_FiscalLimitBinds    = [];
nStates_LafferCurveBinds    = [];
mdl_SetUp                   = {};
arrStates = cell2mat(vertcat(stateRecord{:}));
for iMdl = 1:nMdl
    % Reaches the peak of the Laffer curve
    nStates_LafferCurveBinds(iMdl)  = sum(arrStates(:,iMdl) == 2) + sum(arrStates(:,iMdl) == 4);
    
    % Reaches the Fiscal limit
    nStates_FiscalLimitBinds(iMdl)  = sum(arrStates(:,iMdl) == 3) + sum(arrStates(:,iMdl) == 4);
    
    mdl_SetUp{iMdl} = riseSetUpNames{iMdl} ;
end

% Do the constraints bind?
prob_Laffer         = nStates_LafferCurveBinds' ./ (nSims * simul_periods) .* 100 .* 4 ;
prob_FiscalLimit    = nStates_FiscalLimitBinds' ./ (nSims * simul_periods) .* 100 .* 4 ;
formattedNumbers    = arrayfun(@ (n) sprintf("%1.2f\\%%", n), [round(prob_Laffer, 2), round(prob_FiscalLimit, 2)]);
tbConstraintsBind   = array2table(formattedNumbers, ...
                        'RowNames', mdl_SetUp', ...
                        'VariableNames', {'Peak of the Laffer curve', 'Fiscal limit'});
disp(tbConstraintsBind)

tbPath =  strjoin({pathTables 'Simulation' modelName 'DefR' defProcess ['Table - Simulation - Constraints - All Rules.tex']}, filesep);
edu_Table2Latex(tbConstraintsBind, tbPath)

%% Simulate specific sequence of shocks (INCOMPLETE)

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 1;
simul_regime                            = 1:4;
simul_order                             = 1;
simul_pruned                            = false;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 80;
simul_burn                              = 0;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% Simulate the model
mdl = mdlVector;

[mySims,retcode] = simulate_nonlinear(mdl, ...
    'simul_regime',             simul_regime, ...
    'simul_order',              simul_order, ...
    'simul_pruned',             simul_pruned, ...
    'simul_honor_constraints',  simul_honor_constraints, ...
    'simul_honor_constraints_through_switch', simul_honor_constraints_through_switch, ...
    'simul_anticipate_zero',    simul_anticipate_zero, ...
    'simul_periods',            simul_periods, ...
    'simul_burn',               simul_burn, ...
    'simul_frwrd_back_shoot',   simul_frwrd_back_shoot, ...
    'simul_shock_uncertainty',  simul_shock_uncertainty);

%% Simulate the model

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlSimVector = mdlVector;
for iMdl = 1:length(mdlSimVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
        
        %paramsStructTemp.phiC = 1000;
        %paramsStructTemp.rhoA = 0 ;
        %paramsStructTemp.rhoGG = 0 ;
        %paramsStructTemp.rhoGY = 0 ;
        %paramsStructTemp.rhoM = 0 ;
        %paramsStructTemp.sigmaM = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        paramsStructTemp.sigmaM    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlSimVector(iMdl)      = set(mdlSimVector(iMdl),'parameters', paramsStructTemp);
    mdlSimVector(iMdl)      = solve(mdlSimVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 4;
simul_regime                            = 1:4;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 300000;
simul_burn                              = 30000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false
%simul_sig                               = 0; % 0 = return to the deterministic steady state

% Simulate the model
mdl = mdlSimVector;
tic;
simulateModel;
toc;

%%%%%%%%%%%%%% Check convergence
for iMdl = 1:nMdl
    disp(['Model ' num2str(iMdl)]);
    simVals = vertcat(stateRecord{:});
    simVals = cell2mat(simVals(:,iMdl));
    stateChains = reshape(simVals, [simul_periods nSims]);
    
    minReg = min(simVals);
    maxReg = max(simVals);
    
    tConvergence = array2table(NaN(maxReg,nSims));
    for iReg = 1:maxReg
        tConvergence{iReg, :} = sum(stateChains == iReg);
    end
    tConvergence.Properties.RowNames = strcat("Regime ", num2str((1:maxReg)'));
    tConvergence.Properties.VariableNames = strcat("Chain ", num2str((1:nSims)'));
    disp(tConvergence);
end
%%%%%%%%%%%%%%%%%

% Verify regime binding
nStates_FiscalLimitBinds    = [];
nStates_LafferCurveBinds    = [];
mdl_SetUp                   = {};
for iMdl = 1:nMdl
    simVals = vertcat(stateRecord{:});
    simVals = cell2mat(simVals(:,iMdl));
    for iChain = 1:nSims
        % Reaches the peak of the Laffer curve
        nStates_LafferCurveBinds(iMdl)  = sum(simVals == 2) + sum(simVals == 4);

        % Reaches the Fiscal limit
        nStates_FiscalLimitBinds(iMdl)  = sum(simVals == 3) + sum(simVals == 4);
    end
    
    mdl_SetUp{iMdl} = ['Mdl = ', num2str(iMdl), ' GovBondsRiskFree = ' , char(string(mdlVector(iMdl).user_data.conf_govBondsAreRiskFree)), ' Rule = ', mdlVector(iMdl).user_data.conf_policyRule, ', Debt = ', mdlVector(iMdl).user_data.conf_debtLevel];
end

% Do the constraints bind?
prob_Laffer         = nStates_LafferCurveBinds' ./ (nSims * simul_periods) .* 100 ;
prob_FiscalLimit    = nStates_FiscalLimitBinds' ./ (nSims * simul_periods) .* 100 ;
formattedNumbers = arrayfun(@ (n) sprintf("%1.2f%%", n), [round(prob_Laffer, 2), round(prob_FiscalLimit, 2)]);
tbConstraintsBind = array2table(formattedNumbers, ...
                        'RowNames', mdl_SetUp', ...
                        'VariableNames', {'Peak of the Laffer curve binds?', 'Fiscal limit binds?'});
disp(tbConstraintsBind)

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

modelName   = riseModelNames{1};
defProcess  = riseConfigVector(1).conf_defRProcessType;

tbPath =  strjoin({pathTables 'Simulation' ['Table - Simulation - Regimes Frequency - Debt Level ' debtLevel ' - All Rules - ' title_polRule_simple ' - '  title_approxPoint '.tex']}, filesep);
edu_Table2Latex(tbConstraintsBind, tbPath);

%% Simulate the model (specific path)

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlSimVector = mdlVector;
for iMdl = 1:length(mdlSimVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaM    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlSimVector(iMdl)      = set(mdlSimVector(iMdl),'parameters', paramsStructTemp);
    mdlSimVector(iMdl)      = solve(mdlSimVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

nSims                                   = 1;
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 10000;
simul_burn                              = 1000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false
simul_historical_data                   = false;
forecast_cond_exo_vars                   = {};

%%%%%%%%%%%% Specify initial conditions
db = initial_conditions(mdlSimVector(1),1,'steady',1);

% Define Regime Path - We fix the economy in regime 1 all time
db.regime = db.epsA;
db.regime(1, 1, 1:(simul_periods+simul_burn)) = 1;

% Define Shock Path - Null in this case
%db.epsA(1, 1, 1:(simul_periods+simul_burn)) = 0; 

% Choose Active Shocks
shocks = {};

simul_historical_data                   = db;
forecast_cond_exo_vars                  = shocks;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate the model
mdl = mdlSimVector;
tic;
simulateModelHistorical;
toc;

%%%%%%%%%%%%%% Check convergence
for iMdl = 1:nMdl
    disp(['Model ' num2str(iMdl)]);
    simVals = vertcat(stateRecord{:});
    simVals = cell2mat(simVals(:,iMdl));
    stateChains = reshape(simVals, [simul_periods nSims]);
    
    minReg = min(simVals);
    maxReg = max(simVals);
    
    tConvergence = array2table(NaN(maxReg,nSims));
    for iReg = minReg:maxReg
        tConvergence{iReg, :} = sum(stateChains == iReg);
    end
    tConvergence.Properties.RowNames = strcat("Regime ", num2str((minReg:maxReg)'));
    tConvergence.Properties.VariableNames = strcat("Chain ", num2str((1:nSims)'));
    disp(tConvergence);
end
%%%%%%%%%%%%%%%%%

% Verify regime binding
nStates_FiscalLimitBinds    = [];
nStates_LafferCurveBinds    = [];
mdl_SetUp                   = {};
for iMdl = 1:nMdl
    simVals = vertcat(stateRecord{:});
    simVals = cell2mat(simVals(:,iMdl));
    for iChain = 1:nSims
        % Reaches the peak of the Laffer curve
        nStates_LafferCurveBinds(iMdl)  = sum(simVals == 2) + sum(simVals == 4);

        % Reaches the Fiscal limit
        nStates_FiscalLimitBinds(iMdl)  = sum(simVals == 3) + sum(simVals == 4);
    end
    
    mdl_SetUp{iMdl} = ['GovBondsRiskFree = ' , char(string(mdlVector(iMdl).user_data.conf_govBondsAreRiskFree)), ' Rule = ', mdlVector(iMdl).user_data.conf_policyRule, ', Debt = ', mdlVector(iMdl).user_data.conf_debtLevel];
end

% Do the constraints bind?
prob_Laffer         = nStates_LafferCurveBinds' ./ (nSims * simul_periods) .* 100 ;
prob_FiscalLimit    = nStates_FiscalLimitBinds' ./ (nSims * simul_periods) .* 100 ;
formattedNumbers = arrayfun(@ (n) sprintf("%1.2f%%", n), [round(prob_Laffer, 2), round(prob_FiscalLimit, 2)]);
tbConstraintsBind = array2table(formattedNumbers, ...
                        'RowNames', mdl_SetUp', ...
                        'VariableNames', {'Peak of the Laffer curve binds?', 'Fiscal limit binds?'});
disp(tbConstraintsBind)

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

modelName   = riseModelNames{1};
defProcess  = riseConfigVector(1).conf_defRProcessType;

tbPath =  strjoin({pathTables 'Simulation' modelName 'DefR' defProcess ['Table - Simulation - Regimes Frequency - Debt Level ' debtLevel ' - All Rules - ' title_polRule_simple ' - '  title_approxPoint '.tex']}, filesep);
edu_Table2Latex(tbConstraintsBind, tbPath);

%% Simulate models using the same shocks

%%%%%%%%%%%%%% SET-UP MODELS
mdlSimVector = mdlVector;
for iMdl = 1:length(mdlSimVector)
    paramsStructTemp = struct();
    if polRule_simple
        paramsStructTemp.piiBar  = 0 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlSimVector(iMdl)      = set(mdlSimVector(iMdl),'parameters', paramsStructTemp);
    mdlSimVector(iMdl)      = solve(mdlSimVector(iMdl));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Initiate randomizer
rng(1900)

iMdl                                    = 1;
nSims                                   = 1;
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = false;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = 500;
simul_burn                              = 1000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

activeShocks = {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'};
activeShocks = activeShocks(ismember(activeShocks, shock_list{iMdl}));

% Simulate
simRecord       = struct();
stateRecord     = [];
for iSim = 1:nSims
        
        disp(strcat('Simulation ', string(iSim)));
        
        % Run the simulation for the 1st model
        [mySims,simStates,retcode] = simulate(mdlSimVector(1), ...
            'simul_regime',             simul_regime, ...
            'simul_order',              simul_order, ...
            'simul_pruned',             simul_pruned, ...
            'simul_honor_constraints',  simul_honor_constraints, ...
            'simul_honor_constraints_through_switch', simul_honor_constraints_through_switch, ...
            'simul_anticipate_zero',    simul_anticipate_zero, ...
            'simul_periods',            simul_periods, ...
            'simul_burn',               simul_burn, ...
            'simul_frwrd_back_shoot',   simul_frwrd_back_shoot, ...
            'simul_shock_uncertainty',  simul_shock_uncertainty);
        
        % Transform single vector into cell for allowing for N models
        if ~iscell(simStates)
            simStates = {simStates};
        end
        
        % Recover the shocks from the first model
        db = initial_conditions(mdlVector(1), 1, 'steady', 1);
        for iShock = 1:length(shock_list)
            if sum(ismember(shock_list{1}{iShock}, activeShocks)) == 0
                db.(shock_list{1}{iShock})(1, 1, 1:(simul_periods)) = 0;
            else
                db.(shock_list{1}{iShock})(1, 1, 1:(simul_periods)) = ...
                    mySims.(shock_list{1}{iShock}).data(:,1);
            end
        end
        
        % set initial conditions for regime
        % create the field db.regime; by copying db.e_r to it it will have the right
        % character (ts); if the first step is omitted, db.regime will have the
        % character double and the next simulation does not run
        %db.regime = db.e_r; 
        %db.regime(1, 1, 1:301)=mysims.regime.data;
        
        % Redo the simulation
        [mySims,simStates,retcode] = simulate(mdlSimVector, ...
            'simul_regime',             simul_regime, ...
            'simul_order',              simul_order, ...
            'simul_pruned',             simul_pruned, ...
            'simul_honor_constraints',  simul_honor_constraints, ...
            'simul_honor_constraints_through_switch', simul_honor_constraints_through_switch, ...
            'simul_anticipate_zero',    simul_anticipate_zero, ...
            'simul_periods',            simul_periods, ...
            'simul_burn',               simul_burn, ...
            'simul_frwrd_back_shoot',   simul_frwrd_back_shoot, ...
            'simul_shock_uncertainty',  simul_shock_uncertainty, ...
            'simul_historical_data',    db, ...
            'forecast_cond_exo_vars',   activeShocks);
        
        % Transform single vector into cell for allowing for N models
        if ~iscell(simStates)
            simStates = {simStates};
        end
        
        if iSim == 1
            simRecord = mySims;
            stateRecord = simStates;
        else
            simRecord = [simRecord, mySims];
            stateRecord = [stateRecord; simStates];
        end
        
end


% Verify regime binding
nStates_FiscalLimitBinds    = [];
nStates_LafferCurveBinds    = [];
mdl_SetUp                   = {};
for iMdl = 1:nMdl
    % Reaches the peak of the Laffer curve
    nStates_LafferCurveBinds(iMdl)  = sum(vec([stateRecord{:}{:,iMdl}]) == 2) + sum(vec([stateRecord{:}{:,iMdl}]) == 4);
    
    % Reaches the Fiscal limit
    nStates_FiscalLimitBinds(iMdl)  = sum(vec([stateRecord{:}{:,iMdl}]) == 3) + sum(vec([stateRecord{:}{:,iMdl}]) == 4);
    
    mdl_SetUp{iMdl} = ['Policy rule = ', mdlSimVector(iMdl).user_data.conf_policyRule, ', Debt Level = ', mdlSimVector(iMdl).user_data.conf_debtLevel];
end

% Do the constraints bind?
tbConstraintsBind = table(nStates_LafferCurveBinds', nStates_FiscalLimitBinds', ...
                        'RowNames', mdl_SetUp', ...
                        'VariableNames', {'Peak of the Laffer curve Binds?', 'Fiscal Limit Binds?'});
disp(tbConstraintsBind)

%% Save simulations

save([pathSaved filesep 'sims_1' '.mat'], 'simRecord', 'stateRecord');

%% Simulation: PLOT state transitions

close all;

refMdl = 1;
varlist = {'polDef', 'rGap', 'rPolicy', 'Pii', 'nrPolicy', 'y', 'c', 'n'};
plotConfig;
modelName = riseModelNames{refMdl};
defProcess = riseConfigVector(1).conf_defRProcessType;

arrStates = cell2mat(vertcat(stateRecord{:}));
nPlotPeriods_sim = 1:size(arrStates, 1);

f = figure;

nMdl    = size(arrStates, 2);
minState = 0; %min(vec(cell2mat(stateRecord{:})));
maxState = max(vec(arrStates));

% Set layout and reduce empty spaces in the graph
tl = tiledlayout(nMdl, 1);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';

previousInterpreter = edu_GraphSetInterpreter('latex');
for iMdl = 1:nMdl
    tl_i = nexttile();
    yPlot = arrStates(:, iMdl);
    p = stairs(yPlot, '-.ob');
    p = edu_GraphSetInterpreter('latex');
    p.LineStyle     = ':';
    p.Color         = 'blue';
    p.LineWidth     = 2;
    p.Title         = 'States';
    p.XLabel        = 'time';
    p.XLimits       = [nPlotPeriods_sim(1) nPlotPeriods_sim(end)];
    p.Title         = setUpNames{iMdl};
    ylim([minState maxState]);
end

lg  = legend(tl_i, 'Orientation', 'Horizontal', 'NumColumns', nMdl); 
%lg.Position = 'North'; % <-- Legend placement with tiled layout

set(tl.Children, 'FontSize', 16);
set(f, 'Position',  [100, 100, 800, 800]); % resize figure
edu_GraphSetInterpreter(previousInterpreter);

exportgraphics(f, ...
    strjoin({pathImages 'IRFs' modelName 'DefR' defProcess ['Graph - Simulation - States' ' - DebtLevel ' debtLevel '.png']}, filesep));

%% Simulation: PLOT selected variables

close all;

% Periods (regime transition)
mdlToFindRegime     = 1;
regimeToShow        = 2;
regimePeriods       = find(simStates{mdlToFindRegime} == regimeToShow);
if isempty(regimePeriods)
    nPlotPeriods_sim    = 1:simul_periods;
else
    nPlotPeriods_sim    = (regimePeriods(1)-50):(regimePeriods(1)+50);
end

refMdl = 1;
plotConfig;
modelName = riseModelNames{refMdl};
defProcess = riseConfigVector(1).conf_defRProcessType;

f = figure;

% Set layout and reduce empty spaces in the graph
tl = tiledlayout(nLins, nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';

previousInterpreter = edu_GraphSetInterpreter('latex');
for iPlot = 1:length(varlist)
    tl_i = nexttile();

    yPlot = mySims.(varlist{iPlot}).values;
    p = stackedplot(yPlot);
    
    p = edu_GraphSetInterpreterStackedPlot(p, 'latex');
    
    
    p.LineStyle     = '-';
    p.Color         = 'blue';
    p.LineWidth     = 2;
    p.Title         = vtexNames{iPlot};
    p.XLabel        = 'time';
    p.XLimits       = [nPlotPeriods_sim(1) nPlotPeriods_sim(end)];
    
    if mod(iPlot, nCols) == 1
        p.DisplayLabels = setUpNames;
    else
        p.DisplayLabels = repmat({''},1, size(yPlot,2));
    end
    
end
set(tl.Children, 'FontSize', 16);
set(f, 'Position',  [100, 100, 800, 800]); % resize figure
edu_GraphSetInterpreter(previousInterpreter);

exportgraphics(f, ...
    strjoin({pathImages 'IRFs' modelName 'DefR' defProcess ['Graph - Simulation - All Variables' ' - DebtLevel ' debtLevel '.png']}, filesep));

%% Simulation: histogram of selected variables (unconditional distribution)

%%%%%%%%%%%
graphType = 'Rules Comparison'; % 'Rules Comparison', 'Risk-Free vs Risky', 'Welfare 2nd Order'
%%%%%%%%%%%

histMdl = mdlSimVector;
histVars = {'nrPolicy', 'Pii', 'bY', 'r1fisLim_1_2', 'y', 'c', 'cR', 'cNR', 'n', 'gY', 'tau', 'tax'};
histTexVars = {'$i_t$', '$\Pi_t$', '$\frac{B_t}{Y_t}$', '$\mathcal{D}_t$', '$Y_t$', '$C_t$', '$C^R_t$', '$C^{NR}_t$', '$N_t$', '$\frac{G_t}{Y_t}$', '$\tau_t$', '$T_t$'};

N = 4;
C = linspecer(N);
plotColors = {C(1,:), C(2,:), C(3,:), C(4,:)};
plotStyles = {'-', ':', '--', '-.'};
plotWidths = {1, 2, 1, 1};
nLins = 3;
nCols  = ceil(length(histVars)/nLins);

fontSize = 12;

f = figure;
set(gca,'LooseInset',get(gca,'TightInset'));
previousInterpreter = edu_GraphSetInterpreter('latex');
tl = tiledlayout(nLins, nCols);
tl.Padding          = 'compact';
tl.TileSpacing      = 'compact';

for iVar = 1:length(histVars)
    nexttile();
    h = gobjects(length(histMdl), 1);
    for iMdl = 1:length(histMdl)
        y_tau = [];
        for iChain = 1:length(simRecord)
            y_tau = [y_tau; simRecord(iChain).(histVars{iVar}).values];
        end
        
        % Rescale values
        if strcmp(histVars{iVar}, 'nrPolicy')
            histTexVars{iVar} = '$i_t$ (\% annualized)';
            y_tau(:,iMdl) = ((1 + y_tau(:,iMdl)).^4 - 1) .* 100;
        elseif strcmp(histVars{iVar}, 'Pii')
            histTexVars{iVar} = '$\pi_t$ (\% annualized)';
            y_tau(:,iMdl) = (y_tau(:,iMdl).^4 - 1) .* 100;
        elseif strcmp(histVars{iVar}, 'bY')
            histTexVars{iVar} = '$\frac{B_t}{4 Y_t}$ (\%)';
            y_tau(:,iMdl) = (y_tau(:,iMdl)./4) .* 100;
        elseif strcmp(histVars{iVar}, 'fisLim_1_2')
            histTexVars{iVar} = '$\mathcal{D}_t$ (\%)';
            y_tau(:,iMdl) = y_tau(:,iMdl) .* 100;
        elseif strcmp(histVars{iVar}, 'tau')
            histTexVars{iVar} = '$\tau_t$ (\%)';
            y_tau(:,iMdl) = y_tau(:,iMdl) .* 100;
        elseif strcmp(histVars{iVar}, 'gY')
            histTexVars{iVar} = '$\frac{G_t}{Y_t}$ (\%)';
            y_tau(:,iMdl) = y_tau(:,iMdl) .* 100;
        end
        
        % Check whether there is variance
        if std(y_tau(:,iMdl)) > 0
            dist = fitdist(y_tau(:,iMdl),'Kernel','Kernel','epanechnikov');
            %if strcmp(histVars{iVar}, 'nrPolicy')
            %    xRange = icdf(dist, 0.01):((max(y(:,iMdl)) - icdf(dist, 0.0001))/100):icdf(dist, 0.999);
            %else
                xRange = min(y_tau(:,iMdl)):((max(y_tau(:,iMdl)) - min(y_tau(:,iMdl)))/100):max(y_tau(:,iMdl)); 
            %end
            hold on;
            h(iMdl) = area([xRange], pdf(dist, xRange), 'LineStyle', plotStyles{iMdl}, 'DisplayName', riseSetUpNames{iMdl});
            h(iMdl).LineWidth = plotWidths{iMdl};
            h(iMdl).FaceAlpha = 0.2;
            
            %h(iMdl).Color = plotColors{iMdl};
            xl = xline(median(y_tau(:,iMdl)), ['--']);
            xl.Color =  h(iMdl).FaceColor;
            xl.LineWidth = 1.5;
        else
            hold on;
            h(iMdl) = xline(median(y_tau(:,iMdl)), ['-']);
            h(iMdl).Color = plotColors{iMdl};
            h(iMdl).LineWidth = 1.5;
        end
        set(gca, 'FontSize', fontSize);
    end
    title(histTexVars{iVar});
    axis tight;
    if iVar == 1
       lg = legend(h, 'Location', 'best', 'FontSize', fontSize, 'Orientation', 'horizontal');
       lg.Layout.Tile = 'North';
    end
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

title_graphType = graphType;
if strcmp(graphType, 'Risk-Free vs Risky')
    title_debtLevel = 'Both';
else
    if strcmp(graphType, 'Welfare 2nd Order')
       title_graphType = [title_graphType ' ' histMdl(1).user_data.conf_policyRule ' and ' histMdl(2).user_data.conf_policyRule];
    end
    
    if strcmp(debtLevel, 'High')
        title_debtLevel = ['High ' num2str(histMdl(1).user_data.conf_defPolTarget*100)];
    else
        title_debtLevel = debtLevel;
    end
end


set(f, 'Position',  [100, 0, 1200, nLins*250]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Distribution' title_stickyPrices graphType ['Graph - Distribution - ' title_graphType ' - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Simulation: histogram of nrPolicy and pii variables (unconditional distribution)

%%%%%%%%%%%%%%%%%%
removeNegativeNrPolicy = false;
graphType = 'Rules Comparison'; % 'Rules Comparison' or 'Risk-Free vs Risky'
plotVerticalBars = [1,2];
%%%%%%%%%%%%%%%%%%

histMdl = mdlSimVector;
histVars = {'nrPolicy', 'Pii'}; % 'nrPolicy', 'Pii', 'rPolicy', 'rGov'
histTexVars = {'$i_t$', '$\Pi_t$'}; % '$i_t$', '$\Pi_t$', '$r_t$', '$r^{Gov}_t$'
plotColors = {'b', 'r', 'k', [0.4660 0.6740 0.1880]};
plotStyles = {'-', ':', '-.', '--'};
plotLineWidth = {1, 3, 1, 1};
plotXLims = {[-5 45], [0 10]}; % [-5 45], [0 10]
nLins = 1;
nCols  = ceil(length(histVars)/nLins);

fontSize = 12;
tol = 1e-12;

f = figure;
previousInterpreter = edu_GraphSetInterpreter('latex');
tl = tiledlayout(nLins, nCols);
tl.Padding          = 'compact';
tl.TileSpacing      = 'compact';

yNrPolicy   = [];
yPii        = [];
yRPolicy    = [];
yRGov       = [];
yRNa       = [];
yREF       = [];
for iChain = 1:length(simRecord)
    yNrPolicy = [yNrPolicy; simRecord(iChain).('nrPolicy').values];
    yPii = [yPii; simRecord(iChain).('Pii').values];
    yRPolicy = [yRPolicy; simRecord(iChain).('rPolicy').values];
    yRGov = [yRGov; simRecord(iChain).('rGov').values];
    yRNa = [yRNa; simRecord(iChain).('rNa').values];
    yREF = [yREF; simRecord(iChain).('rEF').values];
end

for iVar = 1:length(histVars)
    nexttile();
    h = gobjects(length(histMdl), 1);
    for iMdl = 1:length(histMdl)
        
        
        % Remove negative interest rates
        if removeNegativeNrPolicy
            yPii_Mdl        = yPii(yNrPolicy(:,iMdl) >= 0, iMdl);
            yNrPolicy_Mdl   = yNrPolicy(yNrPolicy(:,iMdl) >= 0, iMdl);
            yRPolicy_Mdl    = yRPolicy(yNrPolicy(:,iMdl) >= 0, iMdl);
            yRGov_Mdl       = yRGov(yNrPolicy(:,iMdl) >= 0, iMdl);
            yRNa_Mdl        = yRNa(yNrPolicy(:,iMdl) >= 0, iMdl);
            yREF_Mdl        = yREF(yNrPolicy(:,iMdl) >= 0, iMdl);
        else
            yPii_Mdl        = yPii(:, iMdl);
            yNrPolicy_Mdl   = yNrPolicy(:, iMdl);
            yRPolicy_Mdl    = yRPolicy(:, iMdl);
            yRGov_Mdl       = yRGov(:, iMdl);
            yRNa_Mdl        = yRNa(:, iMdl);
            yREF_Mdl        = yREF(:, iMdl);
        end
        
        if strcmp(histVars{iVar}, 'nrPolicy')
            histTexVars{iVar} = '$i_t$ (\% annualized)';
            y_tau = yNrPolicy_Mdl;
        elseif strcmp(histVars{iVar}, 'Pii')
            histTexVars{iVar} = '$\pi_t$ (\% annualized)';
            y_tau = yPii_Mdl;
        elseif strcmp(histVars{iVar}, 'rPolicy')
            histTexVars{iVar} = '$r_t$ (\% annualized)';
            y_tau = yRPolicy_Mdl;
        elseif strcmp(histVars{iVar}, 'rGov')
            histTexVars{iVar} = '$r^{Gov}_t$ (\% annualized)';
            y_tau = yRGov_Mdl;
        elseif strcmp(histVars{iVar}, 'rNa')
            histTexVars{iVar} = '$r^{n}_t$ (\% annualized)';
            y_tau = yRNa_Mdl;
        elseif strcmp(histVars{iVar}, 'rEF')
            histTexVars{iVar} = '$r^{Eff}_t$ (\% annualized)';
            y_tau = yREF_Mdl;
        end
        
        if ismember(histVars{iVar}, {'nrPolicy', 'rPolicy', 'rGov', 'rNa', 'rEF'})
            yData = ((1 + y_tau).^4 - 1).*100;
        elseif strcmp(histVars{iVar}, 'Pii')
            yData = ((y_tau).^4 - 1).*100;
        end
        
        % Check whether there is variance
        if std(yData) > tol
            
            dist = fitdist(yData,'Kernel','Kernel','epanechnikov');
            %xRange = icdf(dist, 0.1):((max(yData) - icdf(dist, 0.0001))/100):icdf(dist, 0.999);
            xRange = min(yData):((max(yData) - min(yData))/100):max(yData);
            hold on;
            h(iMdl) = plot([xRange], pdf(dist, xRange), 'LineStyle', plotStyles{iMdl}, 'DisplayName', riseSetUpNames{iMdl});
            h(iMdl).Color = plotColors{iMdl};
            h(iMdl).LineWidth = plotLineWidth{iMdl};
            
            if ismember(iMdl, plotVerticalBars)
                xl = xline(mean(yData), ['--' plotColors{iMdl}]);
                xl.LineWidth = 2;
                
                yLim = get(gca,'ylim');
                xLim = get(gca,'xlim');
                if iMdl == plotVerticalBars(1)
                    txt = ['$\leftarrow$ ' num2str(mean(yData), '%.1f')];
                    text(mean(yData), yLim(1) + (yLim(2)-yLim(1))/10, txt, 'HorizontalAlignment', 'left', 'FontSize',12)
                elseif iMdl == plotVerticalBars(2)
                    txt = [num2str(mean(yData), '%.1f') '$\rightarrow$ '];
                    text(mean(yData), yLim(1) + (yLim(2)-yLim(1))/10, txt, 'HorizontalAlignment', 'right', 'FontSize',12)
                end
            end
        else
            hold on;
            %h(iMdl) = xline(median(yData), ['-' plotColors{iMdl}]);
        end
    end
   
    set(gca, 'FontSize', fontSize);
    title(histTexVars{iVar});
    xlim(plotXLims{iVar});
    
    %axis tight;
    if iVar == 1
       lg = legend(h, 'Location', 'best', 'Orientation', 'horizontal', 'FontSize', fontSize);
       lg.Layout.Tile = 'North';
       ylabel('probability density');
    end
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

if strcmp(graphType, 'Risk-Free vs Risky')
    title_debtLevel = 'Both';
else
    if strcmp(debtLevel, 'High')
        title_debtLevel = ['High ' num2str(histMdl(1).user_data.conf_defPolTarget*100)];
    else
        title_debtLevel = debtLevel;
    end
end

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if removeNegativeNrPolicy
   title_ZLB = 'ZLB'; 
else
    title_ZLB = 'No ZLB'; 
end

set(f, 'Position',  [100, 0, 1200, nLins*500]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Distribution' title_stickyPrices graphType ['Graph - Distribution - ' histVars{1} ' and ' histVars{2} ' - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint ' - ' title_ZLB '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);


%%%%%% Compare DSS with the ergodic mean
ssMatrix     = [];
ssVarNames   = [];
for iMdl = 1:nMdl
    nRegimes     = histMdl(iMdl).markov_chains.regimes_number;
    ssMatrix     = [ssMatrix, [histMdl(iMdl).solution.ss{:}]];
    ssVarNames   = [ssVarNames; strcat("Mdl ", num2str(iMdl), "; Regime ", string(1:nRegimes)')];
end
tSS = array2table(ssMatrix);
tSS.Properties.VariableNames    = ssVarNames;
tSS.Properties.RowNames         = mdlVector(1).endogenous.name;
disp(tSS);

nMdl = length(histMdl);
nReg = histMdl(1).markov_chains.regimes_number;

tComp = array2table(NaN(nMdl, 6));
tComp.Properties.VariableNames = {'DSS', 'EM', 'Bias', 'DSS ', 'EM ', 'Bias '};
tComp.Properties.RowNames = riseSetUpNames;

tComp2 = array2table(NaN(nMdl, 6));
tComp2.Properties.VariableNames = {'DSS', 'EM', 'Bias', 'DSS ', 'EM ', 'Bias '};
tComp2.Properties.RowNames = riseSetUpNames;

tComp3 = array2table(NaN(nMdl, 6));
tComp3.Properties.VariableNames = {'DSS', 'EM', 'Bias', 'DSS ', 'EM ', 'Bias '};
tComp3.Properties.RowNames = riseSetUpNames;

for iMdl = 1:length(histMdl)
    iNrPolicy   = find(strcmp(tSS.Properties.RowNames, 'nrPolicy'));
    iPii        = find(strcmp(tSS.Properties.RowNames, 'Pii'));
    iRPolicy    = find(strcmp(tSS.Properties.RowNames, 'rPolicy'));
    iRGov       = find(strcmp(tSS.Properties.RowNames, 'rGov'));
    iRNa        = find(strcmp(tSS.Properties.RowNames, 'rNa'));
    iREF        = find(strcmp(tSS.Properties.RowNames, 'rEF'));
    
    %%% TABLE 1: nrPolicy and Pii
    tComp{iMdl,1} = ((tSS.(['Mdl ', num2str(iMdl), '; Regime 1'])(iNrPolicy) + 1)^4 - 1)*100;
    tComp{iMdl,1} = round(tComp{iMdl,1}, 1);
    
    if removeNegativeNrPolicy
        tComp{iMdl,2} = mean((((1 + yNrPolicy(yNrPolicy(:,iMdl) >= 0, iMdl)).^4) - 1).*100);
        tComp{iMdl,5} = mean(((yPii(yNrPolicy(:,iMdl) >= 0, iMdl).^4) - 1).*100);
    else
        tComp{iMdl,2} = mean((((1 + yNrPolicy(:,iMdl)).^4) - 1).*100);
        tComp{iMdl,5} = mean(((yPii(:,iMdl).^4) - 1).*100);
    end
    tComp{iMdl,2} = round(tComp{iMdl,2}, 1);
    tComp{iMdl,5} = round(tComp{iMdl,5}, 1);
    
    tComp{iMdl,3} = tComp{iMdl,2} - tComp{iMdl,1};
    
    tComp{iMdl,4} = (tSS.(['Mdl ', num2str(iMdl), '; Regime 1'])(iPii)^4 - 1)*100;
    tComp{iMdl,4} = round(tComp{iMdl,4}, 1);
    
    tComp{iMdl,6} = tComp{iMdl,5} - tComp{iMdl,4};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% TABLE 2: rPolicy and rGov
    tComp2{iMdl,1} = ((tSS.(['Mdl ', num2str(iMdl), '; Regime 1'])(iRPolicy) + 1)^4 - 1)*100;
    tComp2{iMdl,1} = round(tComp2{iMdl,1}, 1);
    
    if removeNegativeNrPolicy
        tComp2{iMdl,2} = mean((((1 + yRPolicy(yNrPolicy(:,iMdl) >= 0, iMdl)).^4) - 1).*100);
        tComp2{iMdl,5} = mean((((1 + yRGov(yNrPolicy(:,iMdl) >= 0, iMdl)).^4) - 1).*100);
    else
        tComp2{iMdl,2} = mean((((1 + yRPolicy(:,iMdl)).^4) - 1).*100);
        tComp2{iMdl,5} = mean((((1 + yRGov(:,iMdl)).^4) - 1).*100);
    end
    tComp2{iMdl,2} = round(tComp2{iMdl,2}, 1);
    tComp2{iMdl,5} = round(tComp2{iMdl,5}, 1);
    
    tComp2{iMdl,3} = tComp2{iMdl,2} - tComp2{iMdl,1};
    
    tComp2{iMdl,4} = ((tSS.(['Mdl ', num2str(iMdl), '; Regime 1'])(iRGov) + 1)^4 - 1)*100;
    tComp2{iMdl,4} = round(tComp2{iMdl,4}, 1);
    
    tComp2{iMdl,6} = tComp2{iMdl,5} - tComp2{iMdl,4};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% TABLE 3: rNa and rEF
    tComp3{iMdl,1} = ((tSS.(['Mdl ', num2str(iMdl), '; Regime 1'])(iRNa) + 1)^4 - 1)*100;
    tComp3{iMdl,1} = round(tComp3{iMdl,1}, 1);
    
    if removeNegativeNrPolicy
        tComp3{iMdl,2} = mean((((1 + yRNa(yNrPolicy(:,iMdl) >= 0, iMdl)).^4) - 1).*100);
        tComp3{iMdl,5} = mean((((1 + yREF(yNrPolicy(:,iMdl) >= 0, iMdl)).^4) - 1).*100);
    else
        tComp3{iMdl,2} = mean((((1 + yRNa(:,iMdl)).^4) - 1).*100);
        tComp3{iMdl,5} = mean((((1 + yREF(:,iMdl)).^4) - 1).*100);
    end
    tComp3{iMdl,2} = round(tComp3{iMdl,2}, 1);
    tComp3{iMdl,5} = round(tComp3{iMdl,5}, 1);
    
    tComp3{iMdl,3} = tComp3{iMdl,2} - tComp3{iMdl,1};
    
    tComp3{iMdl,4} = ((tSS.(['Mdl ', num2str(iMdl), '; Regime 1'])(iREF) + 1)^4 - 1)*100;
    tComp3{iMdl,4} = round(tComp3{iMdl,4}, 1);
    
    tComp3{iMdl,6} = tComp3{iMdl,5} - tComp3{iMdl,4};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
disp(tComp);
disp(tComp2);
disp(tComp3);

% In LaTex define a new column type in the beginnig of the document: \newcolumntype{Y}{>{\centering\arraybackslash}X}
colAlignment = '|YYY|YYY';
tabWidth = '1.0\\textwidth';
colNames = '\multicolumn{3}{c}{$i_t$ (\% annualized)} & \multicolumn{3}{|c}{$\pi_t$ (\% annualized)} \\ & DSS & EM & Bias & DSS & EM & Bias';
tbPath =  strjoin({pathTables 'Simulation' ['Table - Distribution - nrPolicy and Pii - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint ' - ' title_ZLB '.tex']}, filesep);
edu_Table2Latex(tComp, tbPath, 'colAlignment', colAlignment, 'colNames', colNames, 'tabWidth', tabWidth);

% In LaTex define a new column type in the beginnig of the document: \newcolumntype{Y}{>{\centering\arraybackslash}X}
colAlignment = '|YYY|YYY';
tabWidth = '1.0\\textwidth';
colNames = '\multicolumn{3}{c}{$r_t$ (\% annualized)} & \multicolumn{3}{|c}{$r^{Gov}_t$ (\% annualized)} \\ & DSS & EM & Bias & DSS & EM & Bias';
tbPath =  strjoin({pathTables 'Simulation' ['Table - Distribution - rPolicy and rGov - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint ' - ' title_ZLB '.tex']}, filesep);
edu_Table2Latex(tComp2, tbPath, 'colAlignment', colAlignment, 'colNames', colNames, 'tabWidth', tabWidth);

% In LaTex define a new column type in the beginnig of the document: \newcolumntype{Y}{>{\centering\arraybackslash}X}
colAlignment = '|YYY|YYY';
tabWidth = '1.0\\textwidth';
colNames = '\multicolumn{3}{c}{$r^n_t$ (\% annualized)} & \multicolumn{3}{|c}{$r^{Eff}_t$ (\% annualized)} \\ & DSS & EM & Bias & DSS & EM & Bias';
tbPath =  strjoin({pathTables 'Simulation' ['Table - Distribution - rNa and rEF - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint ' - ' title_ZLB '.tex']}, filesep);
edu_Table2Latex(tComp3, tbPath, 'colAlignment', colAlignment, 'colNames', colNames, 'tabWidth', tabWidth);

%% Simulation:Stylized histogram of Pii

%%%%%%%%%%%%%%%%%%
removeNegativeNrPolicy = false;
graphType = 'Stylized';
plotVerticalBars = [1,2,3,4];
%%%%%%%%%%%%%%%%%%

histMdl = mdlSimVector;
histVars = {'Pii'}; % 'nrPolicy', 'Pii', 'rPolicy', 'rGov'
histTexVars = {'$\Pi_t$'}; % '$i_t$', '$\Pi_t$', '$r_t$', '$r^{Gov}_t$'
plotColors = {'b', 'r', 'k', [0.4660 0.6740 0.1880]};
plotStyles = {'-', ':', '-.', '--'};
plotLineWidth = {1, 3, 1, 1};
plotXLims = {[3 6]}; % [-5 45], [0 10]
nLins = 1;
nCols  = ceil(length(histMdl)/nLins);

fontSize = 12;
tol = 1e-9;

f = figure;
previousInterpreter = edu_GraphSetInterpreter('latex');
tl = tiledlayout(nLins, nCols);
tl.Padding          = 'compact';
tl.TileSpacing      = 'compact';

yNrPolicy   = [];
yPii        = [];
yRPolicy    = [];
yRGov       = [];
yRNa       = [];
yREF       = [];
for iChain = 1:length(simRecord)
    yNrPolicy = [yNrPolicy; simRecord(iChain).('nrPolicy').values];
    yPii = [yPii; simRecord(iChain).('Pii').values];
    yRPolicy = [yRPolicy; simRecord(iChain).('rPolicy').values];
    yRGov = [yRGov; simRecord(iChain).('rGov').values];
    yRNa = [yRNa; simRecord(iChain).('rNa').values];
    yREF = [yREF; simRecord(iChain).('rEF').values];
end

for iVar = 1:length(histVars)
    
    h = gobjects(length(histMdl), 1);
    for iMdl = 1:length(histMdl)
        nexttile();
        
        % Remove negative interest rates
        if removeNegativeNrPolicy
            yPii_Mdl        = yPii(yNrPolicy(:,iMdl) >= 0, iMdl);
            yNrPolicy_Mdl   = yNrPolicy(yNrPolicy(:,iMdl) >= 0, iMdl);
            yRPolicy_Mdl    = yRPolicy(yNrPolicy(:,iMdl) >= 0, iMdl);
            yRGov_Mdl       = yRGov(yNrPolicy(:,iMdl) >= 0, iMdl);
            yRNa_Mdl        = yRNa(yNrPolicy(:,iMdl) >= 0, iMdl);
            yREF_Mdl        = yREF(yNrPolicy(:,iMdl) >= 0, iMdl);
        else
            yPii_Mdl        = yPii(:, iMdl);
            yNrPolicy_Mdl   = yNrPolicy(:, iMdl);
            yRPolicy_Mdl    = yRPolicy(:, iMdl);
            yRGov_Mdl       = yRGov(:, iMdl);
            yRNa_Mdl        = yRNa(:, iMdl);
            yREF_Mdl        = yREF(:, iMdl);
        end
        
        if strcmp(histVars{iVar}, 'nrPolicy')
            histTexVars{iVar} = '$i_t$ (\% annualized)';
            y_tau = yNrPolicy_Mdl;
        elseif strcmp(histVars{iVar}, 'Pii')
            histTexVars{iVar} = '$\pi_t$ (\% annualized)';
            y_tau = yPii_Mdl;
        elseif strcmp(histVars{iVar}, 'rPolicy')
            histTexVars{iVar} = '$r_t$ (\% annualized)';
            y_tau = yRPolicy_Mdl;
        elseif strcmp(histVars{iVar}, 'rGov')
            histTexVars{iVar} = '$r^{Gov}_t$ (\% annualized)';
            y_tau = yRGov_Mdl;
        elseif strcmp(histVars{iVar}, 'rNa')
            histTexVars{iVar} = '$r^{n}_t$ (\% annualized)';
            y_tau = yRNa_Mdl;
        elseif strcmp(histVars{iVar}, 'rEF')
            histTexVars{iVar} = '$r^{Eff}_t$ (\% annualized)';
            y_tau = yREF_Mdl;
        end
        
        if ismember(histVars{iVar}, {'nrPolicy', 'rPolicy', 'rGov', 'rNa', 'rEF'})
            yData = ((1 + y_tau).^4 - 1).*100;
        elseif strcmp(histVars{iVar}, 'Pii')
            yData = ((y_tau).^4 - 1).*100;
        end
        
        % Check whether there is variance
        if std(yData) > tol
            
            dist = fitdist(yData,'Kernel','Kernel','epanechnikov');
            %xRange = icdf(dist, 0.1):((max(yData) - icdf(dist, 0.0001))/100):icdf(dist, 0.999);
            xRange = min(yData):((max(yData) - min(yData))/100):max(yData);
            hold on;
            h(iMdl) = plot([xRange], pdf(dist, xRange), 'LineStyle', plotStyles{iMdl}, 'DisplayName', riseSetUpNames{iMdl});
            h(iMdl).Color = plotColors{iMdl};
            h(iMdl).LineWidth = plotLineWidth{iMdl};
            
        else
            hold on;
            
            %h(iMdl) = xline(median(yData), ['-' plotColors{iMdl}]);
        end
    
        if ismember(iMdl, plotVerticalBars)
            xl = xline(mean(yData), [plotStyles{iMdl} plotColors{iMdl}], 'DisplayName', riseSetUpNames{iMdl});
            xl.LineWidth = 2;

            yLim = get(gca,'ylim');
            xLim = get(gca,'xlim');
            if iMdl == plotVerticalBars(1)
                txt = ['$\leftarrow$ ' num2str(mean(yData), '%.1f') ' = $\overline{\pi}$'];
                text(mean(yData), yLim(1) + (yLim(2)-yLim(1))/10, txt, 'HorizontalAlignment', 'left', 'FontSize',12)
            elseif iMdl == plotVerticalBars(2)
                txt = [num2str(mean(yData), '%.1f') '$\rightarrow$ '];
                text(mean(yData), yLim(1) + (yLim(2)-yLim(1))/10, txt, 'HorizontalAlignment', 'right', 'FontSize',12)
            elseif iMdl == plotVerticalBars(3)
                txt = ['$\leftarrow$ ' num2str(mean(yData), '%.1f')];
                text(mean(yData), yLim(1) + (yLim(2)-yLim(1))/10, txt, 'HorizontalAlignment', 'right', 'FontSize',12)
            elseif iMdl == plotVerticalBars(4)
                txt = ['$\leftarrow$ ' num2str(mean(yData), '%.1f') ' = $\overline{\pi}$'];
                text(mean(yData), yLim(1) + (yLim(2)-yLim(1))/10, txt, 'HorizontalAlignment', 'left', 'FontSize',12)
            end

            h(iMdl) = xl;
        end
        
        set(gca, 'FontSize', fontSize);
        %title(histTexVars{iVar});
        xlabel(histTexVars{iVar});
        xlim(plotXLims{iVar});
        
        if iMdl == 1
           ylabel('probability density'); 
        end
        
    end
    
    %axis tight;
    if iVar == 1
       lg = legend(h, 'Location', 'best', 'Orientation', 'horizontal', 'FontSize', fontSize, 'box', 'off');
       lg.Layout.Tile = 'North';
    end
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

if strcmp(graphType, 'Risk-Free vs Risky')
    title_debtLevel = 'Both';
else
    if strcmp(debtLevel, 'High')
        title_debtLevel = ['High ' num2str(histMdl(1).user_data.conf_defPolTarget*100)];
    else
        title_debtLevel = debtLevel;
    end
end

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if removeNegativeNrPolicy
   title_ZLB = 'ZLB'; 
else
    title_ZLB = 'No ZLB'; 
end

set(f, 'Position',  [100, 0, 1200, 300]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Distribution' title_stickyPrices graphType ['Graph - Stylized Dist - ' histVars{1} ' - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint ' - ' title_ZLB '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Simulation: histogram of welfare variables (unconditional distribution)

graphType = 'Welfare Distribution';

histMdl = mdlSimVector;
histVars = {'welfare', 'welfareR', 'welfareNR'};
histTexVars = {'$U_t$', '$U^{R}_t$', '$U^{NR}_t$'};

N = 4;
C = linspecer(N);
plotColors = {C(1,:), C(2,:), C(3,:), C(4,:)};
plotStyles = {'-', ':', '--', '-.'};
plotWidths = {1, 2, 1, 1};
plotXLims = {[-230 -130], [-230 -130], [-230 -130]}; % [-5 45], [0 10]

nLins = 1;
nCols  = ceil(length(histVars)/nLins);

fontSize = 12;

f = figure;
previousInterpreter = edu_GraphSetInterpreter('latex');
tl = tiledlayout(nLins, nCols);
tl.Padding          = 'compact';
tl.TileSpacing      = 'compact';

for iVar = 1:length(histVars)
    nexttile();
    h = gobjects(length(histMdl), 1);
    for iMdl = 1:length(histMdl)
        y_tau = [];
        for iChain = 1:length(simRecord)
            y_tau = [y_tau; simRecord(iChain).(histVars{iVar}).values];
        end
        
        % Check whether there is variance
        if std(y_tau(:,iMdl)) > 0
            dist = fitdist(y_tau(:,iMdl),'Kernel','Kernel','epanechnikov');
            xRange = min(y_tau(:,iMdl)):((max(y_tau(:,iMdl)) - min(y_tau(:,iMdl)))/100):max(y_tau(:,iMdl));
            hold on;
            h(iMdl) = area([xRange], pdf(dist, xRange), 'LineStyle', plotStyles{iMdl}, 'DisplayName', riseSetUpNames{iMdl});
            h(iMdl).LineWidth = plotWidths{iMdl};
            h(iMdl).FaceAlpha = 0.2;
            %h(iMdl).Color = plotColors{iMdl};
            xl = xline(median(y_tau(:,iMdl)), ['--']);
            xl.Color =  h(iMdl).FaceColor;
            xl.LineWidth = 1.5;
        else
            hold on;
            h(iMdl) = xline(median(y_tau(:,iMdl)), ['-']);
            h(iMdl).Color = plotColors{iMdl};
            h(iMdl).LineWidth = 1.5;
        end
        set(gca, 'FontSize', fontSize);
    end
    title(histTexVars{iVar});
    %axis tight;
    xlim(plotXLims{iVar});
    
    if iVar == 1
       lg = legend(h, 'Location', 'best', 'Orientation', 'horizontal', 'FontSize', fontSize);
       lg.Layout.Tile = 'North';
    end
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

if strcmp(debtLevel, 'High')
    title_debtLevel = ['High ' num2str(histMdl(1).user_data.conf_defPolTarget*100)];
else
    title_debtLevel = debtLevel;
end

title_graphType = graphType;
if strcmp(graphType, 'Welfare Distribution')
   title_graphType = [title_graphType ' ' histMdl(1).user_data.conf_policyRule];
end

if strcmp(debtLevel, 'High')
    title_debtLevel = ['High ' num2str(histMdl(1).user_data.conf_defPolTarget*100)];
else
    title_debtLevel = debtLevel;
end

set(f, 'Position',  [100, 0, 1200, nLins*250]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Welfare' title_stickyPrices ['Graph - ' title_graphType ' - ' title_polRule_simple ' - ' title_approxPoint '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Simulation: histogram of selected variables (regime-specific distribution)

histMdl = mdlSimVector;
histVars = {'nrPolicy', 'Pii', 'bY', 'fisLim_1_2', 'y', 'c', 'cR', 'cNR', 'n', 'gY', 'tau', 'tax', 'welfare', 'welfareR', 'welfareNR'};
histTexVars = {'$i_t$', '$\Pi_t$', '$\frac{B_t}{Y_t}$', '$\mathcal{D}_t$', '$Y_t$', '$C_t$', '$C^R_t$', '$C^{NR}_t$', '$N_t$', '$\frac{G_t}{Y_t}$', '$\tau_t$', '$T_t$', '$U_t$', '$U^{R}_t$', '$U^{NR}_t$'};
plotColors = {'b', 'r', 'k', [0.4660 0.6740 0.1880]};
plotStyles = {'-', ':', '-.', '--'};
nLins = 3;
nCols  = ceil(length(histVars)/nLins);

fontSize = 12;

tol = 10e-10;

for iReg = 1:2

    f = figure;
    previousInterpreter = edu_GraphSetInterpreter('latex');
    tl = tiledlayout(nLins, nCols);
    tl.Padding          = 'compact';
    tl.TileSpacing      = 'compact';

    for iVar = 1:length(histVars)
        nexttile();
        h = gobjects(length(histMdl), 1);
        for iMdl = 1:length(histMdl)
            y_tau = [];
            yReg = [];
            for iChain = 1:length(simRecord)
                y_tau       = [y_tau; simRecord(iChain).(histVars{iVar}).values];
                yReg    = [yReg; simRecord(iChain).('regime').values];
            end
            
            yRegIdx = yReg(:,iMdl) == iReg;
            
            % Check whether there is variance
            if std(y_tau(yRegIdx,iMdl)) > tol
                dist = fitdist(y_tau(yRegIdx,iMdl),'Kernel','Kernel','epanechnikov');
                xRange = min(y_tau(yRegIdx,iMdl)):((max(y_tau(yRegIdx,iMdl)) - min(y_tau(yRegIdx,iMdl)))/100):max(y_tau(yRegIdx,iMdl));
                hold on;
                h(iMdl) = plot([xRange], pdf(dist, xRange), 'LineStyle', plotStyles{iMdl}, 'DisplayName', riseSetUpNames{iMdl});
                h(iMdl).Color = plotColors{iMdl};
                xl = xline(median(y_tau(yRegIdx,iMdl)), ['--' plotColors{iMdl}]);
                xl.LineWidth = 1.5;
            elseif std(y_tau(yRegIdx,iMdl)) <= tol
                hold on;
                h(iMdl) = xline(median(y_tau(yRegIdx,iMdl)), ['-' plotColors{iMdl}]);
                h(iMdl).LineWidth = 1.5;
            else
                % No observations of this regime
            end
            set(gca, 'FontSize', fontSize);
        end
        title(histTexVars{iVar});
        axis tight;
        if iVar == 1
           lg = legend(h, 'Location', 'best', 'FontSize', fontSize);
           lg.Layout.Tile = 'North';
        end
    end

    if polRule_simple
       title_polRule_simple = 'Simplified Rule';
    else
       title_polRule_simple = 'Full Rule';
    end

    if steady_state_unique
        title_approxPoint = 'Ergodic Mean';
    else
        title_approxPoint = 'Regime-Specific';
    end

    set(f, 'Position',  [100, 0, 1200, nLins*250]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'Distribution' ['Graph - Distribution Regime ' num2str(iReg) ' - Selected Variables - Debt Level ' debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint '.png']}, filesep));
    edu_GraphSetInterpreter(previousInterpreter);

end

%% Simulation: histogram of transition probabilities

histMdl = mdlSimVector;
histVars = {'taxLim_1_2', 'taxLim_2_1', 'fisLim_1_2', 'fisLim_2_1'};
histTexVars = {'$\text{TL}12_t$', '$\text{TL}21_t$', '$\text{FL}12_t$', '$\text{FL}21_t$'};
plotColors = {'b', 'r', 'k'};

nLins = 1;
nCols  = ceil(length(histVars)/nLins);

f = figure;
previousInterpreter = edu_GraphSetInterpreter('latex');
tl = tiledlayout(nLins, nCols);
%tl.Padding          = 'compact';
tl.TileSpacing      = 'compact';

for iVar = 1:length(histVars)
    nexttile();
    h = gobjects(length(histMdl), 1);
    for iMdl = 1:length(histMdl)
        y_tau = [];
        for iChain = 1:length(simRecord)
            y_tau = [y_tau; simRecord(iChain).(histVars{iVar}).values];
        end
        
        if std(y_tau(yRegIdx,iMdl)) > tol
            dist = fitdist(y_tau(:,iMdl),'Kernel','Kernel','epanechnikov');
            xRange = min(y_tau(:,iMdl)):((max(y_tau(:,iMdl)) - min(y_tau(:,iMdl)))/100):max(y_tau(:,iMdl));
            hold on;
            h = plot([xRange], pdf(dist, xRange));
            h.Color = plotColors{iMdl};
            xl = xline(mean(y_tau(:,iMdl)), [':' plotColors{iMdl}]);
            xl.LineWidth = 1.5;
        else
            hold on;
            h(iMdl) = xline(median(y_tau(:,iMdl)), ['-' plotColors{iMdl}]);
            h(iMdl).LineWidth = 1.5;
        end
        
    end
    title(histTexVars{iVar});
    axis tight;
    currYLim = ylim; 
    ylim([0 min(1, currYLim(2))]);
    if iVar == 1
       legend(riseSetUpNames, 'Location', 'best'); 
    end
end

set(f, 'Position',  [100, 0, 1400, nLins*400]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Distribution' ['Graph - Distribution - Transition Probabilities.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Simulation: correlation of selected variables (unconditional distribution)

histMdl = mdlSimVector;
histVars = {'Pii', 'fisLim_1_2'};
histTexVars = {'$\Pi_t$', '$\mathcal{D}_{t+1}$'};
plotColors = {'b', 'r', 'k'};
plotStyles = {'-', ':', '-.'};
nLins = 1;
nCols  = length(mdlVector);

fontSize = 12;

f = figure;
previousInterpreter = edu_GraphSetInterpreter('latex');
tl = tiledlayout(nLins, nCols);
tl.Padding          = 'compact';
tl.TileSpacing      = 'compact';

y1 = [];
yReg = [];
for iChain = 1:length(simRecord)
    y1 = [y1; simRecord(iChain).(histVars{1}).values];
    yReg = [yReg; simRecord(iChain).('regime').values];
end

y2 = [];
for iChain = 1:length(simRecord)
    y2 = [y2; simRecord(iChain).(histVars{2}).values];
end

yRegIdx = 1:length(yReg(:,iMdl));

% Iterate models
h = gobjects(length(mdlVector), 1);
for iMdl = 1:length(histMdl)
    nexttile();
    
    % Rescale values
    iVar = 1;
    if strcmp(histVars{iVar}, 'Pii')
        histTexVars{iVar} = '$\pi_t$ (\% annualized)';
        y1(:,iMdl) = (y1(:,iMdl).^4 - 1) .* 100;
    end
    iVar = 2;
    if strcmp(histVars{iVar}, 'fisLim_1_2')
        histTexVars{iVar} = '$\mathcal{D}_t$ (\%)';
        y2(:,iMdl) = y2(:,iMdl) .* 100;
    end
    
    % Select sample to plot
    sampleSize = 10000;
    selIdx = randsample(length(y1(yRegIdx,iMdl)), sampleSize);
    
    % Check whether there is variance
    y1_sample = y1(yRegIdx,iMdl);
    y2_sample = y2(yRegIdx,iMdl);
    h(iMdl) = scatter(y1_sample(selIdx), y2_sample(selIdx));
    h(iMdl).Marker = '.';
    h(iMdl).MarkerFaceAlpha = 0.1;
    
    corrVal = corr(y1(yRegIdx,iMdl), y2(yRegIdx,iMdl));
    ylabel(histTexVars{2});
    xlabel({histTexVars{1}, ['Correlation: ' num2str(corrVal, '%.2f')]});
    
    ls = lsline();
    ls.Color = 'black';

    %h(iMdl).Color = plotColors{iMdl};
    set(gca, 'FontSize', fontSize);
    title({riseSetUpNames{iMdl};''});
end
linkaxes(tl.Children);

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

if steady_state_unique
    title_approxPoint = 'Ergodic Mean';
else
    title_approxPoint = 'Regime-Specific';
end

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

if strcmp(debtLevel, 'High')
    title_debtLevel = ['High ' num2str(histMdl(1).user_data.conf_defPolTarget*100)];
else
    title_debtLevel = debtLevel;
end

set(f, 'Position',  [100, 0, 1200, nLins*350]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Correlation' title_stickyPrices ['Graph - Correlation All Regimes - Selected Variables - Debt Level ' title_debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Simulation: correlation of selected variables (regime-specific distribution)

histMdl = mdlSimVector;
histVars = {'nrPolicy', 'fisLim_1_2'};
histTexVars = {'$i_t$', '$\Pi_t$'};
plotColors = {'b', 'r', 'k'};
plotStyles = {'-', ':', '-.'};
nLins = 1;
nCols  = length(mdlVector);

fontSize = 12;

for iReg = 1:4

    f = figure;
    previousInterpreter = edu_GraphSetInterpreter('latex');
    tl = tiledlayout(nLins, nCols);
    tl.Padding          = 'compact';
    tl.TileSpacing      = 'compact';

    y1 = [];
    yReg = [];
    for iChain = 1:length(simRecord)
        y1 = [y1; simRecord(iChain).(histVars{1}).values];
        yReg = [yReg; simRecord(iChain).('regime').values];
    end

    y2 = [];
    for iChain = 1:length(simRecord)
        y2 = [y2; simRecord(iChain).(histVars{2}).values];
    end

    yRegIdx = yReg(:,iMdl) == iReg;
        
    h = gobjects(length(mdlVector), 1);
    for iMdl = 1:length(histMdl)
        nexttile();
        % Check whether there is variance
        h(iMdl) = scatter(y1(yRegIdx,iMdl), y2(yRegIdx,iMdl));
        h(iMdl).Marker = '.';
        h(iMdl).MarkerFaceAlpha = 0.2;
        ls = lsline();
        ls.Color = 'black';
        %h(iMdl).Color = plotColors{iMdl};
        set(gca, 'FontSize', fontSize);
        title(riseSetUpNames{iMdl});
    end
    linkaxes(tl.Children);

    if polRule_simple
       title_polRule_simple = 'Simplified Rule';
    else
       title_polRule_simple = 'Full Rule';
    end

    if steady_state_unique
        title_approxPoint = 'Ergodic Mean';
    else
        title_approxPoint = 'Regime-Specific';
    end

    set(f, 'Position',  [100, 0, 1000, nLins*350]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'Correlation' ['Graph - Correlation Regime ' num2str(iReg) ' - Selected Variables - Debt Level ' debtLevel ' - ' title_polRule_simple ' - ' title_approxPoint '.png']}, filesep));
    edu_GraphSetInterpreter(previousInterpreter);

end

%% Simulation: PLOT tau_t vs. B_t

%%%%% SPECIFY TIME RANGE TO PLOT
rngPeriods = [735 775]; % use a wide range to find the regime-swicthing periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlimRng = [1 (rngPeriods(2)-rngPeriods(1))];

mdl     = 1;
nLins   = 2;
nCols   = 1;

N = 5;
C = linspecer(N);
plotColors = {C(1,:), C(2,:), C(3,:), C(4,:), C(5,:)};
plotStyles = {'-', ':', '--', '-.'};
plotWidths = {2, 3, 2, 2};
fontSize = 14;

f = figure;
set(gca,'LooseInset',get(gca,'TightInset'));
previousInterpreter = edu_GraphSetInterpreter('latex');
% Set layout and reduce empty spaces in the graph
tl = tiledlayout(nLins, nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';

for iMdl = mdl
    nexttile();
    
    y_y         = [];
    y_aTilde    = [];
    y_ddelta     = [];
    y_tau       = [];
    y_tauSdw    = [];
    y_tax       = [];
    y_b         = [];
    y_bY        = [];
    y_nrPolicy  = [];
    y_rPolicy   = [];
    y_nrGov     = [];
    y_iota      = [];
    y_Pii       = [];
    
    y_epsA      = [];
    y_epsG      = [];
    y_epsM      = [];
    
    y_regime    = [];
    y_reg_1_1   = [];
    y_reg_1_2   = [];
    y_reg_1_3   = [];
    y_reg_1_4   = [];
    y_reg_2_1   = [];
    y_reg_2_2   = [];
    y_reg_2_3   = [];
    y_reg_2_4   = [];
    y_reg_3_1   = [];
    y_reg_3_2   = [];
    y_reg_3_3   = [];
    y_reg_3_4   = [];
    y_reg_4_1   = [];
    y_reg_4_2   = [];
    y_reg_4_3   = [];
    y_reg_4_4   = [];
    y_r1fisLim    = [];
    y_r2taxLim    = [];
    for iChain = 1:length(simRecord)
        y_y        = [y_y; simRecord(iChain).('y').data(:,iMdl)];
        y_aTilde   = [y_aTilde; simRecord(iChain).('aTilde').data(:,iMdl)];
        y_ddelta   = [y_ddelta; simRecord(iChain).('ddelta').data(:,iMdl)];
        y_tau      = [y_tau; simRecord(iChain).('tau').data(:,iMdl)];
        y_tauSdw   = [y_tauSdw; simRecord(iChain).('tauSdw').data(:,iMdl)];
        y_tax      = [y_tax; simRecord(iChain).('tax').data(:,iMdl)];
        y_b        = [y_b; simRecord(iChain).('b').data(:,iMdl)];
        y_bY       = [y_bY; simRecord(iChain).('bY').data(:,iMdl)];
        y_nrPolicy = [y_nrPolicy; simRecord(iChain).('nrPolicy').data(:,iMdl)];
        y_rPolicy  = [y_rPolicy; simRecord(iChain).('rPolicy').data(:,iMdl)];
        y_nrGov    = [y_nrGov; simRecord(iChain).('nrGov').data(:,iMdl)];
        y_iota     = [y_iota; simRecord(iChain).('iota').data(:,iMdl)];
        y_Pii      = [y_Pii; simRecord(iChain).('Pii').data(:,iMdl)];
        
        y_epsA = [y_epsA; simRecord(iChain).('epsA').data(:,iMdl)];
        y_epsG = [y_epsG; simRecord(iChain).('epsG').data(:,iMdl)];
        y_epsM = [y_epsM; simRecord(iChain).('epsM').data(:,iMdl)];
        
        y_regime = [y_regime; simRecord(iChain).('regime').data(:,iMdl)];
        
        y_reg_1_1 = [y_reg_1_1; simRecord(iChain).('regime_1_1').data(:,iMdl)];
        y_reg_1_2 = [y_reg_1_2; simRecord(iChain).('regime_1_2').data(:,iMdl)];
        y_reg_1_3 = [y_reg_1_3; simRecord(iChain).('regime_1_3').data(:,iMdl)];
        y_reg_1_4 = [y_reg_1_4; simRecord(iChain).('regime_1_4').data(:,iMdl)];
        
        y_reg_2_1 = [y_reg_2_1; simRecord(iChain).('regime_2_1').data(:,iMdl)];
        y_reg_2_2 = [y_reg_2_2; simRecord(iChain).('regime_2_2').data(:,iMdl)];
        y_reg_2_3 = [y_reg_2_3; simRecord(iChain).('regime_2_3').data(:,iMdl)];
        y_reg_2_4 = [y_reg_2_4; simRecord(iChain).('regime_2_4').data(:,iMdl)];
        
        y_reg_3_1 = [y_reg_3_1; simRecord(iChain).('regime_3_1').data(:,iMdl)];
        y_reg_3_2 = [y_reg_3_2; simRecord(iChain).('regime_3_2').data(:,iMdl)];
        y_reg_3_3 = [y_reg_3_3; simRecord(iChain).('regime_3_3').data(:,iMdl)];
        y_reg_3_4 = [y_reg_3_4; simRecord(iChain).('regime_3_4').data(:,iMdl)];
        
        y_reg_4_1 = [y_reg_4_1; simRecord(iChain).('regime_4_1').data(:,iMdl)];
        y_reg_4_2 = [y_reg_4_2; simRecord(iChain).('regime_4_2').data(:,iMdl)];
        y_reg_4_3 = [y_reg_4_3; simRecord(iChain).('regime_4_3').data(:,iMdl)];
        y_reg_4_4 = [y_reg_4_4; simRecord(iChain).('regime_4_4').data(:,iMdl)];
        
        y_r1fisLim = [y_r1fisLim; simRecord(iChain).('r1fisLim').data(:,iMdl)];
        y_r2taxLim = [y_r2taxLim; simRecord(iChain).('r2taxLim').data(:,iMdl)];
    end
    y_tauMax = get(mdlVector(iMdl),'parameters').tauMaxBar(1) .* ones(length(y_tau),1);
    
    %%%%% Transform variables
    y_bY = y_bY / 4 * 100;
    y_tau = y_tau * 100;
    y_tauSdw = y_tauSdw * 100;
    y_tauMax = y_tauMax * 100;
    %%%%%%%%%%%%%%%%
    
    yPlot = [y_tau, y_tauSdw, y_tauMax];
    p = plot(yPlot(rngPeriods(1):rngPeriods(2), :));
    p(1).LineStyle = plotStyles{1};
    p(1).Color = plotColors{1};
    p(1).LineWidth = plotWidths{1};
    p(2).LineStyle = plotStyles{2};
    p(2).Color = plotColors{2};
    p(2).LineWidth = plotWidths{2};
    p(3).LineStyle = plotStyles{3};
    p(3).Color = 'black';
    p(3).LineWidth = plotWidths{3};
    ylabel('$\%$');
    ylim([30 55]);
    hold on;
    yyaxis right;
    p2 = plot(y_bY(rngPeriods(1):rngPeriods(2), :));
    p2.LineStyle = plotStyles{4};
    p2.Color = plotColors{4};
    p2.LineWidth = plotWidths{4};
    lg = legend('$\tau_t$', '$\tau^{sdw}_t$', '$\tau^{max}_t$', '$\frac{B_t}{4Y_t} \,\%$ (right axis)', ...
        'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', fontSize);
    lg.Layout.Tile = 'North';
    xlabel('time');
    ylabel('$\%$');
    xlim(xlimRng);
    %ylim([  -0.02 + min(vec(yPlot)), ... 
    %        +0.02 + max(vec(yPlot))]);
    set(gca, 'FontSize', fontSize);
    
%     nexttile();
%     yPlot = [y_tau, y_bY];
%     %yPlot = movavg(yPlot,'linear',20);
%     yyaxis left;
%     p = plot(yPlot(:,1));
%     p.LineStyle = '-';
%     p.Color = 'blue';
%     p.LineWidth = 2;
%     ylabel('$\tau_t \,\%$');
%     yyaxis right;
%     hold on;
%     p = plot(yPlot(:,2));
%     p.LineStyle = ':';
%     p.Color = 'green';
%     p.LineWidth = 2;
%     ylabel('$\frac{B_t}{4Y_t} \,\%$');
%     xlabel('time');
%     xlim(rngPeriods);
%     %ylim([  -0.02 + min(vec(yPlot)), ... 
%     %        +0.02 + max(vec(yPlot))]);
%     set(gca, 'FontSize', 12);

    nexttile();
    
    y_peakProb  = (y_regime == 1) .* y_reg_1_2 + ...
                  (y_regime == 2) .* y_reg_2_2 + ...
                  (y_regime == 3) .* y_reg_3_2 + ...
                  (y_regime == 4) .* y_reg_4_2 + ...
                  (y_regime == 1) .* y_reg_1_4 + ...
                  (y_regime == 2) .* y_reg_2_4 + ...
                  (y_regime == 3) .* y_reg_3_4 + ...
                  (y_regime == 4) .* y_reg_4_4;
    
    y_defProb   = (y_regime == 1) .* y_reg_1_3 + ...
                  (y_regime == 2) .* y_reg_2_3 + ...
                  (y_regime == 3) .* y_reg_3_3 + ...
                  (y_regime == 4) .* y_reg_4_3 + ...
                  (y_regime == 1) .* y_reg_1_4 + ...
                  (y_regime == 2) .* y_reg_2_4 + ...
                  (y_regime == 3) .* y_reg_3_4 + ...
                  (y_regime == 4) .* y_reg_4_4;
              
    p1 = plot(100 .* y_defProb(rngPeriods(1):rngPeriods(2), :));
    p1.Color = 'black';
    p1.LineWidth = 2;
    %hold on;
    %a = area( 100.* (y_fisLim == 1) .* (y_taxLim == 1) );
    %a.FaceAlpha = 0.2;
    hold on;
    y = 100.* (y_regime == 2);
    a = bar( y(rngPeriods(1):rngPeriods(2), :), 'BarWidth', 1);
    a.FaceAlpha = 0.2;
    a.EdgeColor = 'none';
    hold on;
    y = 100.* (y_regime == 3);
    a = bar( y(rngPeriods(1):rngPeriods(2), :), 'BarWidth', 1);
    a.FaceAlpha = 0.2;
    a.EdgeColor = 'none';
    hold on;
    y = 100.* (y_regime == 4);
    a = bar( y(rngPeriods(1):rngPeriods(2), :), 'BarWidth', 1);
    a.FaceAlpha = 0.2;
    a.EdgeColor = 'none';
    
    lg = legend('Prob. of default binding', ...
                'Regime 2', 'Regime 3', 'Regime 4', ...
        'Location', 'northoutside', ...
        'Orientation', 'horizontal', 'FontSize', fontSize);
    lg.Layout.Tile = 'south';
    xlabel('time');
    ylabel('$\%$');
    xlim(xlimRng);
    ylim([0, 100*max(y_defProb(rngPeriods(1):rngPeriods(2), :))]);
    set(gca, 'FontSize', fontSize);
    
end

set(f, 'Position',  [100, 300, 800, 600]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Laffer Curve' ['Graph - Simulated Laffer Curve.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% SIMULATION MOMENTS

varlist = {'y', 'nrPolicy', 'Pii'};

% Calculate moments
simMean     = struct();
simVariance = struct();
for iVar = 1:length(varlist)
    simMean.(varlist{iVar})     = mean([simRecord.(varlist{iVar})]);
    simVariance.(varlist{iVar}) = var([simRecord.(varlist{iVar})]);
end

% Capture different configs
tblVars         = {'y', 'nrPolicy', 'Pii'};
tblVarsNames    = {'$Y_t$', '$i_t$', '$\Pi$'};
meanData    = [];
varData     = [];
for iMdl = 1:nMdl
    for iVar = 1:length(tblVars)
        %idx = iMdl:nMdl:(length(tblVars)*nMdl);
        meanData    = [meanData; mean(simMean.(tblVars{iVar})(iMdl))];
        varData     = [varData; mean(simVariance.(tblVars{iVar})(iMdl))];
    end
end

% Display results
meanDataStr     = arrayfun(@(x) num2str(x,2), meanData, 'UniformOutput', false);
varDataStr      = arrayfun(@(x) num2str(x,2), varData, 'UniformOutput', false);
tCombs          = cell2table([repmat(tblVarsNames',nMdl,1), meanDataStr, varDataStr]);
tCombs.Properties.VariableNames = {'Variable', 'Mean', 'Variance'};
tCombs.Properties.VariableNames = {'Variable', 'Mean', 'Variance'};
disp(tCombs);

%save('welfare.mat', 'meanData', 'varData', 'simMean', 'simVariance');

%%%%%%%%%%%%%%%%%% CHECK
idxStart    = max(1, find(y_regime == 2, 1) - 20);
idxEnd      = min(length(y_regime), find(y_regime == 2, 1) + 20);
disp(['Regime is binding at period ' num2str(find(y_regime == 2, 1))]);

tData = [       (idxStart:idxEnd)', ...
                y_regime(idxStart:idxEnd), ...
                y_tauSdw(idxStart:idxEnd), ...
                y_tau(idxStart:idxEnd), ...
                y_tax(idxStart:idxEnd), ...
                y_b(idxStart:idxEnd), ...
                y_aTilde(idxStart:idxEnd), ...
                y_y(idxStart:idxEnd), ...
                y_bY(idxStart:idxEnd), ...
                y_nrPolicy(idxStart:idxEnd), ...
                y_rPolicy(idxStart:idxEnd), ...
                y_nrGov(idxStart:idxEnd), ...
                y_Pii(idxStart:idxEnd), ...
                y_epsA(idxStart:idxEnd), ...
                y_epsG(idxStart:idxEnd), ...
                y_epsM(idxStart:idxEnd), ...
        ];

tSim = array2table(tData);
tSim.Properties.VariableNames = {'Time', 'Regime', '$\tau^{sdw}$', '$\tau$', '$T_t$', '$B_t$', '$A_t$', '$Y_t$', '$\frac{B_t}{Y_t}$', '$i_t$', '$r_t$', 'i^{Gov}_t', '\Pi_t', '$\varepsilon^A_t$', '$\varepsilon^G_t$', '$\varepsilon^M_t$'};
disp(tSim)

%% Simulation: PLOT simulated nrPolicy vs. Pii

mdl     = mdlVector;
nLins   = 2;
nCols   = 2;

f = figure;
oldInterpreter = edu_GraphSetInterpreter('latex');
% Set layout and reduce empty spaces in the graph
tl = tiledlayout(nLins, nCols);
%tl.TileSpacing  = 'compact';
%tl.Padding      = 'compact';

previousInterpreter = edu_GraphSetInterpreter('latex');
for iMdl = 1:length(mdl)
    nexttile();
    
    y_nrPolicy  = [];
    y_Pii       = [];
    for iChain = 1:length(simRecord)
        y_nrPolicy  = [y_nrPolicy; simRecord(iChain).('nrPolicy').data(:,iMdl)];
        y_Pii       = [y_Pii; simRecord(iChain).('Pii').data(:,iMdl)];
    end
    %y_tauMax = get(mdlVector(iMdl),'parameters').tauMaxBar(1) .* ones(length(y_tau),1);
    
    y_Plot_nrPolicy  = randsample(y_nrPolicy, 10000);
    y_Plot_Pii       = randsample(y_Pii, 10000);
    p = scatter(y_Plot_Pii, y_Plot_nrPolicy, 1);
    %p.LineWidth = 0.1;
    %p.MarkerEdgeColor = 'b';
    %p.MarkerFaceColor = [0 0.5 0.5];
    %p(1).LineStyle = '-';
    %p(1).Color = 'blue';
    %p(1).LineWidth = 2;
    %p(2).LineStyle = '--';
    %p(2).Color = 'red';
    %p(2).LineWidth = 2;
    %legend('$\tau_t$', '$\tau^{exp}_t$', '$\tau^{max}_t$', 'Location', 'best');
    xlabel('$\Pi_t$');
    ylabel('$i_t$');
    %ylim([  -0.02 + min(vec(yPlot)), ... 
    %        +0.02 + max(vec(yPlot))]);
    
    hold on;
    p = polyfit(y_Pii, y_nrPolicy, 1);
    x1 = linspace(min(y_Pii), max(y_Pii), 100);
    y1 = polyval(p, x1);    
    plot(x1, y1, 'LineWidth', 3);
    
    hold on;
    plot(mean(y_Pii), mean(y_nrPolicy), 'k.', 'markersize', 20)
    
    title(riseSetUpNames{iMdl});
    set(gca, 'FontSize', 12);
end
linkaxes(tl.Children);

edu_GraphSetInterpreter(previousInterpreter);
set(f, 'Position',  [100, 300, 800, 600]); % resize figure

%% Simulation: PLOT simulated Laffer Curve

mdl     = mdlVector;
nLins   = 2;
nCols   = 2;

f = figure;
oldInterpreter = edu_GraphSetInterpreter('latex');
% Set layout and reduce empty spaces in the graph
tl = tiledlayout(nLins, nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';

previousInterpreter = edu_GraphSetInterpreter('latex');
for iMdl = 1:length(mdl)
    nexttile();
    
    y_tau       = [];
    y_tauMax    = [];
    y_tax       = [];
    y_y         = [];
    for iChain = 1:length(simRecord)
        y_y    = [y_tau; simRecord(iChain).('y').data(:,iMdl)];
        y_tau    = [y_tau; simRecord(iChain).('tau').data(:,iMdl)];
        y_tax = [y_tax; simRecord(iChain).('tax').data(:,iMdl)];
    end
    y_tauMax = get(mdlVector(iMdl),'parameters').tauMaxBar(1) .* ones(length(y_tau),1);
    
    y_Plot_tau  = randsample(y_tau, 10000);
    y_Plot_tax  = randsample(y_tax, 10000);
    p = scatter(y_Plot_tau, y_Plot_tax);
    p.LineWidth = 0.1;
    p.MarkerEdgeColor = 'b';
    p.MarkerFaceColor = [0 0.5 0.5];
    %p(1).LineStyle = '-';
    %p(1).Color = 'blue';
    %p(1).LineWidth = 2;
    %p(2).LineStyle = '--';
    %p(2).Color = 'red';
    %p(2).LineWidth = 2;
    %legend('$\tau_t$', '$\tau^{exp}_t$', '$\tau^{max}_t$', 'Location', 'best');
    xlabel('$\tau_t$');
    ylabel('$T_t$');
    %ylim([  -0.02 + min(vec(yPlot)), ... 
    %        +0.02 + max(vec(yPlot))]);
    
    hold on;
    p = polyfit(y_tau, y_tax, 2);
    x1 = linspace(0, max(y_tau), 100);
    y1 = polyval(p, x1);    
    plot(x1, y1, 'LineWidth', 3);
    
    title(riseSetUpNames{iMdl});
    set(gca, 'FontSize', 12);
end
linkaxes(tl.Children);

edu_GraphSetInterpreter(previousInterpreter);
set(f, 'Position',  [100, 300, 800, 600]); % resize figure

%% Build specific IRFs

buildSpecificIRFs;

irfs_RealSector;

%% Compute impulse responses for all models simultaneously

% %%% GIRF mechanics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://github.com/jmaih/RISE_toolbox/issues/31
% you shock one path over all simulation periods
% you shock a second path over all but the first simulation period
% you take the difference between the two
% you repeat the exercise N times
% you average the results
%
% The paths are just regular simulations i.e. all shocks are considered.
% The only difference between two paths is in the first period where only the shock of interest is activated.
% The theory of pruning is not developed for regime switching models and in that world, there are different types of IRFs that can be computed.
% Regime-specific IRFs are computed on the assumption that you stay in the same regime over the entire simulation period although there is a probability of switching
%.Generalized IRFs are automatically triggered when the model is nonlinear, but they can also be called when the model has more than one regime, which makes it nonlinear.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate randomizer
rng(1900)

useParallel = false;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 40;
irf_type                                = 'irf'; % 'irf', 'girf'
simul_honor_constraints                 = true;
simul_honor_constraints_through_switch  = true;
irf_regime_specific                     = true;
irf_draws                               = 10000;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_pruned                            = true;

% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%% Compute impulse responses for all models simultaneously (with Confidence Interval)

% Initiate randomizer
rng(1900)

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlVector(1),  'irf_shock_sign', shockSize);

nSamples                                = 100;
nIRFPeriods                             = 80;
irf_type                                = 'girf'; % 'irf', 'girf'
simul_honor_constraints                 = true;
simul_honor_constraints_through_switch  = true;
irf_regime_specific                     = false;
irf_draws                               = 10000;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_pruned                            = true;

tic;
computeIRFs_CI;
toc;

%% Save IRFs

save([pathSaved filesep 'irfs_girfs_2' '.mat'], 'myirfs');

%% Plot the impulse responses

close all;

refMdl = 1;
regimeId = 1;
consistentAxes = false;

% varlist = { 'aTilde', 'recTilde', 'def', 'g', 'mShock', 'bbeta', 'polDef', 'tau', ...
%                         'y', 'cStar', 'cRStar', 'cNRStar', 'w', 'n', ...
%                         'nrPolicy' , 'rRN', 'rPolicy', 'Pii', 'muFisLim', 'bY', 'bd', 'tax'};

varlist = {'polDef', 'rGap', 'rPolicy', 'Pii', 'nrPolicy', 'y', 'c', 'n'};
%varlist = {'polDef', 'y', 'c', 'n', 'bY'};
                    
plotConfig;

modelName   = riseModelNames{1};
defProcess  = riseConfigVector(1).conf_defRProcessType;
if riseConfigVector(1).conf_govBondsAreRiskFree
    riskyPolicy = 'Risk-Free';
else
    riskyPolicy = 'Risky';
end

nPlotPeriods_irf = 1:40; % nPlotPeriods
modelsToHide = {};

shockVars = {'epsA',        'aTilde';
             'epsBeta',     'bbeta';
             'epsD',        'def';
             'epsG',        'g';
             'epsM',        'mShock';
             'epsPolDef',   'polDef';
             'epsR',        'recTilde';
             'epsTau',      'tau'
             };

%selectedShocks = shock_names{refMdl};
selectedShocks = {  'epsA', ...
                    'epsBeta', ...
                    'epsG', ...
                    'epsM', ...
                    'epsPolDef', ...
                    'epsTau'
                };
         
previousInterpreter = edu_GraphSetInterpreter('latex');
for iShock = 1:numel(selectedShocks)
    
    %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);
    
    if size(endo_ss,1) == 1
        endo_ss_aux = [endo_ss{1,regimeId}];
    else
        endo_ss_aux = endo_ss(:,regimeId)';
        endo_ss_aux = [endo_ss_aux{:,:}];
    end
    endo_ss_names = endo_names;
    [irfLevel, irfLevelChg, irfPercChg] = ...
        structIRFsAsContainers(myirfs.(selectedShocks{iShock}), endo_ss_names, endo_ss_aux, regimeId);
    
    graph = findobj( 'Type', 'Figure', 'Name', selectedShocks{iShock} );
    if isempty(graph)
        graph = [];
    end
    
    % Exclude shock variables not being used
    shockPlotVars = shockVars(~strcmp(shockVars(:,1), selectedShocks{iShock}), 2);
    shockPlotVars = shockPlotVars(~strcmp(shockPlotVars, 'polDef')); % Remove some variables from the exclusion list
    shockPlotVarNames = plotVarNames(~ismember(plotVars, shockPlotVars));
    shockPlotVars = setdiff(plotVars, shockPlotVars, 'stable');   
    
    % Plot IRFs
    graph = plotIRFs(graph, shockPlotVars, nPlotPeriods_irf, suptitleStr, legendStr, ...
        irfLevel, irfLevelChg, irfPercChg, nCols, shockPlotVarNames, ...
        legendFirstSubPlot, setFirstPeriodToNaN, customizeLines);
    graph.Name = selectedShocks{iShock};
    
    % Hide specific models
    for iGraph = 1:length(graph.Children)
        if isa(graph.Children(iGraph), 'matlab.graphics.axis.Axes')
            for iMdl = 1:length(modelsToHide)
                if ~isempty(graph.Children(iGraph).Children)
                    if isa(graph.Children(iGraph).Children(modelsToHide{iMdl}), 'matlab.graphics.chart.primitive.Line')
                        graph.Children(iGraph).Children(modelsToHide{iMdl}).Visible = 'off';
                    end
                end
            end
        end
    end
    
    set(graph.Children.Children, 'FontSize', 12);
    set(graph.Children.Children, 'FontWeight', 'bold');
    set(findobj(graph.Children.Children, 'Type', 'Legend'), 'FontSize', 10);
    
    leg = findobj(graph.Children.Children, 'Type', 'Legend');
    leg.Location = 'northeast';
    
    if consistentAxes
         linkaxes(graph.Children.Children); 
    end
    
    set(gcf, 'Position',  [100, 100, 1100, nLins*200]); % resize figure
    exportgraphics(graph, strjoin({pathImages 'IRFs' modelName 'DefR' defProcess riskyPolicy  ['Graph - ' selectedShocks{iShock} ' - DebtLevel ' debtLevel ' - Shock ' num2str(shockSize) '.png']}, filesep));
end
edu_GraphSetInterpreter(previousInterpreter);

%% Plot GRAPH: Risky 3 Rules

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlGraphVector = mdlVector;
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau    = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 40;
irf_type                                = 'irf'; % 'irf', 'girf'
simul_honor_constraints                 = true;
simul_honor_constraints_through_switch  = true;
irf_regime_specific                     = true;
irf_draws                               = 100000;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_pruned                            = true;

% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%% PLOT
close all;

refMdl = 1;
regimeId = 1;
consistentAxes = false;

% varlist = { 'aTilde', 'recTilde', 'def', 'g', 'mShock', 'bbeta', 'polDef', 'tau', ...
%                         'y', 'cStar', 'cRStar', 'cNRStar', 'w', 'n', ...
%                         'nrPolicy' , 'rRN', 'rPolicy', 'Pii', 'muFisLim', 'bY', 'bd', 'tax'};

varlist = {'rGap', 'rRN', 'rPolicy', 'muFisLim', 'stdFisLim', 'Pii', 'nrPolicy', 'y', 'c', 'n' 'b'};
%varlist = {'polDef', 'y', 'c', 'n', 'bY'};
                    
plotConfig;

modelName   = riseModelNames{1};
defProcess  = riseConfigVector(1).conf_defRProcessType;
if riseConfigVector(1).conf_govBondsAreRiskFree
    riskyPolicy = 'Risk-Free';
else
    riskyPolicy = 'Risky';
end

nPlotPeriods_irf = 1:40; % nPlotPeriods
modelsToHide = {};

shockVars = {'epsA',        'aTilde';
             'epsBeta',     'bbeta';
             'epsD',        'def';
             'epsG',        'g';
             'epsM',        'mShock';
             'epsPolDef',   'polDef';
             'epsR',        'recTilde';
             'epsTau',      'tau'
             };

%selectedShocks = shock_names{refMdl};
selectedShocks = {  'epsA', ...
                    ...%'epsBeta', ...
                    'epsG', ...
                    'epsM', ...
                    ...%'epsPolDef', ...
                    ...%'epsTau'
                };
         
previousInterpreter = edu_GraphSetInterpreter('latex');
for iShock = 1:numel(selectedShocks)
    
    %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);
    
    if size(endo_ss,1) == 1
        endo_ss_aux = [endo_ss{1,regimeId}];
    else
        endo_ss_aux = endo_ss(:,regimeId)';
        endo_ss_aux = [endo_ss_aux{:,:}];
    end
    endo_ss_names = endo_names;
    [irfLevel, irfLevelChg, irfPercChg] = ...
        structIRFsAsContainers(myirfs.(selectedShocks{iShock}), endo_ss_names, endo_ss_aux, regimeId);
    
    graph = findobj( 'Type', 'Figure', 'Name', selectedShocks{iShock} );
    if isempty(graph)
        graph = [];
    end
    
    % Exclude shock variables not being used
    shockPlotVars = shockVars(~strcmp(shockVars(:,1), selectedShocks{iShock}), 2);
    shockPlotVars = shockPlotVars(~strcmp(shockPlotVars, 'polDef')); % Remove some variables from the exclusion list
    shockPlotVarNames = plotVarNames(~ismember(plotVars, shockPlotVars));
    shockPlotVars = setdiff(plotVars, shockPlotVars, 'stable');   
    
    % Plot IRFs
    graph = plotIRFs(graph, shockPlotVars, nPlotPeriods_irf, suptitleStr, legendStr, ...
        irfLevel, irfLevelChg, irfPercChg, nCols, shockPlotVarNames, ...
        legendFirstSubPlot, setFirstPeriodToNaN, customizeLines);
    graph.Name = selectedShocks{iShock};
    
    % Hide specific models
    for iGraph = 1:length(graph.Children)
        if isa(graph.Children(iGraph), 'matlab.graphics.axis.Axes')
            for iMdl = 1:length(modelsToHide)
                if ~isempty(graph.Children(iGraph).Children)
                    if isa(graph.Children(iGraph).Children(modelsToHide{iMdl}), 'matlab.graphics.chart.primitive.Line')
                        graph.Children(iGraph).Children(modelsToHide{iMdl}).Visible = 'off';
                    end
                end
            end
        end
    end
    
    set(graph.Children.Children, 'FontSize', 12);
    set(graph.Children.Children, 'FontWeight', 'bold');
    set(findobj(graph.Children.Children, 'Type', 'Legend'), 'FontSize', 10);
    
    leg = findobj(graph.Children.Children, 'Type', 'Legend');
    leg.Location = 'northeast';
    
    if consistentAxes
         linkaxes(graph.Children.Children); 
    end
    
    set(gcf, 'Position',  [100, 100, 1300, nLins*200]); % resize figure
    exportgraphics(graph, strjoin({pathImages 'IRFs' modelName 'DefR' defProcess riskyPolicy  ['3 rules - IRF Type ' irf_type ' -  Shock ' selectedShocks{iShock} ' - DebtLevel ' debtLevel ' - ShockSize ' num2str(shockSize) '.png']}, filesep));
end
edu_GraphSetInterpreter(previousInterpreter);

%% (BROKEN) Plot assets IRFs to different shocks in the same graph - Real Sector

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlGraphVector = mdlVector;
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau     = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 40;
irf_type                                = 'irf'; % 'irf', 'girf'
irf_regime_specific                     = false;
irf_draws                               = 1;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = nIRFPeriods;
simul_burn                              = 0;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%%%%% PLOT
close all;

refMdl = 1;
varlist = {'polDef', 'y', 'c', 'n', 'bY'};
plotConfig;
modelName   = riseModelNames{refMdl};
defProcess  = riseConfigVector(1).conf_defRProcessType;
policyRule  = riseConfigVector(1).conf_policyRule;

shocksToPlotList        = {};
shocksToPlotNamesList   = {};

shocksToPlotList{1}            = {'epsA', 'epsG', 'epsM'}; % {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'}
shocksToPlotNamesList{1}       = {'$\varepsilon^A$', '$\varepsilon^G$', '$\varepsilon^\mathcal{M}$'};
%shocksToPlotList{2}            = {'epsPolDef'}; % {'epsPolDef', 'epsBeta', 'epsTau'}
%shocksToPlotNamesList{2}       = {'$\varepsilon^{\mathcal{D}}$', '$\varepsilon^\beta$', '$\varepsilon^\tau$'};

plotVars_assets         = {'y', 'c', 'n', 'bY', 'nrPolicy', 'rPolicy', 'Pii'};
plotVarsNames_assets    = {'$Y_t$', '$C_t$', '$N_t$', '${\frac{B}{Y}}_t$', '$i_t$', '$r_t$', '$\Pi_t$'};

nCols = 2 + length(plotVars_assets);
nFigs = length(shocksToPlotList);

fontSize = 12;
pColors = linspecer(3); 
%customizeLines = false;
%customizeLines = struct();
customizeLines.custom_LineWidth = {1, 2, 1, 1, 1, 1, 1};
customizeLines.custom_LineStyle = {'-', ':', '--', '-.', '-', ':', '--'};
customizeLines.custom_LineColor = {pColors(1,:), pColors(2,:), pColors(3,:), 'blue', 'blue', 'blue', 'blue'};
customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};

for iFigure = 1:nFigs
    
    shocksToPlot        = shocksToPlotList{iFigure};
    shocksToPlotNames   = shocksToPlotNamesList{iFigure};
    nShocks = length(shocksToPlot);
    
    idxAxtiveShocks         = ismember(shocksToPlot, fieldnames(myirfs));
    shocksToPlot            = shocksToPlot(idxAxtiveShocks);
    shocksToPlotNames       = shocksToPlotNames(idxAxtiveShocks);
    
    f = figure;
    set(gca, 'FontSize', 14);
    previousInterpreter = edu_GraphSetInterpreter('latex');

    % Set layout and reduce empty spaces in the graph
    tl = tiledlayout(numel(shocksToPlot), nCols);
    tl.TileSpacing  = 'compact';
    tl.Padding      = 'compact';

    for iShock = 1:numel(shocksToPlot)

        %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);

        regimeId = 1;
        if size(endo_ss,1) == 1
            endo_ss_aux = [endo_ss{1,1}]; 
        else
            endo_ss_aux = endo_ss(:,regimeId)';
            endo_ss_aux = [endo_ss_aux{:,:}];
        end
        endo_ss_names = endo_names;
        [irfLevel, irfLevelChg, irfPercChg] = ...
            structIRFsAsContainers(myirfs.(shocksToPlot{iShock}), endo_ss_names, endo_ss_aux, regimeId);

        legendSubPlot = 2;
        tl_h = nexttile();
        tl_h = [tl_h, nexttile()];
        tl_h = [tl_h, nexttile()];
        tl_h = [tl_h, nexttile()];
        tl_h = [tl_h, nexttile()];
        tl_h = [tl_h, nexttile()];
        tl_h = [tl_h, nexttile()];
        tl_h = [tl_h, nexttile()];

        graph = plotIRFs(tl_h, plotVars_assets, nPlotPeriods, suptitleStr, legendStr, ...
                irfLevel, irfLevelChg, irfPercChg, nCols, plotVarsNames_assets, ...
                legendSubPlot, setFirstPeriodToNaN, customizeLines);

        % Format subplots
        for iAx = 3:length(graph)
            edu_GraphDrawZeroAxis(graph(iAx));
            box(graph(iAx));
        end
            
        % Link axes
        %linkaxes(graph);
        dataGraph = cell2mat(arrayfun(@(x) [x.Children.YData], graph(legendSubPlot+1:end), 'UniformOutput', false));
        minDataGraph = min(dataGraph);
        maxDataGraph = max(dataGraph);
        diffDataGraph = maxDataGraph - minDataGraph;
        %if diffDataGraph == 0 && minDataGraph == 0
        %    ylim(graph, [minDataGraph - 1, maxDataGraph + 1]);
        %else
        %    ylim(graph, [minDataGraph - 0.2*diffDataGraph, maxDataGraph + 0.2*diffDataGraph]);
        %end

        % Row label
        text(graph(2), mean(xlim(graph(2)))*0.5, mean(ylim(graph(2))), ...
            shocksToPlotNames{iShock}, ...
            'Interpreter', 'latex', ...
            'FontSize', 20);

        % First row
        if iShock == 1
            graph(1).Legend.Location = 'North';
        end

        % Hide after the first row
        if iShock > 1
            title(graph(1:end), '');
            legend(graph(1), 'off');
        end

    %     % Between lines 1 and 3
    %     if iShock >= 1 && iShock <= 3
    %         set(graph(4:6), 'Color', [0.9 0.9 0.9]);
    %     end

        % Hide after column 2
        %ylabel(graph(4:end), '');

        % Hide after column 2
        %yticks(graph(3:end), []);

        % Turning its axis off
        axis(graph(2), 'off');

    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize legend
    ax = nexttile(1);
    ax.Legend.FontSize = 12;

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 16;
    end
    
    if polRule_simple
       title_polRule_simple = 'Simplified Rule';
    else
       title_polRule_simple = 'Full Rule';
    end
    
    
    set(f, 'Position',  [100, 0, 1200, 150*nShocks]); % resize figure
    exportgraphics(f, ...
        strjoin({pathImages 'IRFs' modelName 'DefR' defProcess 'RF vs Risky' ['Graph - Real Sector - All Shocks - ' 'IRF Type ' irf_type ' - PolicyRule ' policyRule ' ' - ' title_polRule_simple ' num2str(iFigure) '.png']}, filesep));
    edu_GraphSetInterpreter(previousInterpreter);
end

%% Plot assets IRFs to different shocks in the same graph - Real Sector (Regime-Specific)

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlGraphVector = mdlVector;
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0 ;
        %paramsStructTemp.phi     = 1.1 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau     = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 40;
irf_type                                = 'irf'; % 'irf', 'girf'
irf_regime_specific                     = true;
irf_draws                               = 1;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = nIRFPeriods;
simul_burn                              = 0;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%%%%% PLOT
close all;

refMdl = 1;
varlist = {'polDef', 'y', 'c', 'n', 'bY'};
plotConfig;
modelName   = riseModelNames{refMdl};
defProcess  = riseConfigVector(1).conf_defRProcessType;
policyRule  = riseConfigVector(1).conf_policyRule;
stickyPrices = conf_stickyPrices;

shocksToPlotList        = {};
shocksToPlotNamesList   = {};

shocksToPlotList{1}            = {'epsA', 'epsG', 'epsM'}; % {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'}
shocksToPlotNamesList{1}       = {'$\varepsilon^A$', '$\varepsilon^G$', '$\varepsilon^\mathcal{M}$'};
%shocksToPlotList{2}            = {'epsPolDef'}; % {'epsPolDef', 'epsBeta', 'epsTau'}
%shocksToPlotNamesList{2}       = {'$\varepsilon^{\mathcal{D}}$', '$\varepsilon^\beta$', '$\varepsilon^\tau$'};

plotVars_assets         = {'y', 'c', 'n', 'bY', 'nrPolicy', 'rPolicy', 'Pii'};
plotVarsNames_assets    = {'$Y_t$', '$C_t$', '$N_t$', '${\frac{B}{Y}}_t$', '$i_t$', '$r_t$', '$\Pi_t$'};

nCols = length(plotVars_assets);
nFigs = length(shocksToPlotList);

fontSize = 12;
pColors = linspecer(3); 
%customizeLines = false;
%customizeLines = struct();
customizeLines.custom_LineWidth = {1, 1, 2, 1, 1, 1, 1};
customizeLines.custom_LineStyle = {'-', '--', ':', '-.', '-', '-', '-'};
customizeLines.custom_LineColor = {'blue', 'red', 'black', 'blue', 'blue', 'blue', 'blue'};
customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};

nReg = nRegimes;
for regimeId = 1:nReg
    for iFigure = 1:nFigs

        shocksToPlot        = shocksToPlotList{iFigure};
        shocksToPlotNames   = shocksToPlotNamesList{iFigure};
        nShocks = length(shocksToPlot);

        idxAxtiveShocks         = ismember(shocksToPlot, fieldnames(myirfs));
        shocksToPlot            = shocksToPlot(idxAxtiveShocks);
        shocksToPlotNames       = shocksToPlotNames(idxAxtiveShocks);

        f = figure;
        set(gca, 'FontSize', 14);
        previousInterpreter = edu_GraphSetInterpreter('latex');

        % Set layout and reduce empty spaces in the graph
        tl = tiledlayout(numel(shocksToPlot), nCols);
        tl.TileSpacing  = 'compact';
        tl.Padding      = 'compact';
        
        for iShock = 1:numel(shocksToPlot)

            %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);

            if size(endo_ss,1) == 1
                endo_ss_aux = [endo_ss{1,1}]; 
            else
                endo_ss_aux = endo_ss(:,regimeId)';
                endo_ss_aux = [endo_ss_aux{:,:}];
            end
            endo_ss_names = endo_names;
            [irfLevel, irfLevelChg, irfPercChg] = ...
                structIRFsAsContainers(myirfs.(shocksToPlot{iShock}), endo_ss_names, endo_ss_aux, regimeId);

            legendSubPlot = 1;
            tl_h = nexttile();
            for iVar = 2:length(plotVars_assets)
                tl_h = [tl_h, nexttile()];
            end

            graph = plotIRFs(tl_h, plotVars_assets, nPlotPeriods, suptitleStr, [], ...
                    irfLevel, irfLevelChg, irfPercChg, nCols, plotVarsNames_assets, ...
                    0, setFirstPeriodToNaN, customizeLines);
            
            ylabel(tl_h(1), {shocksToPlotNames{iShock};tl_h(1).YLabel.String}, 'FontSize', 20, 'Interpreter', 'latex');
                
            % Format subplots
            for iAx = 1:length(graph)
                edu_GraphDrawZeroAxis(graph(iAx));
                box(graph(iAx));
                legend('off');
            end

            % Link axes
            %linkaxes(graph);
            dataGraph = cell2mat(arrayfun(@(x) [x.Children.YData], graph(legendSubPlot:end), 'UniformOutput', false));
            minDataGraph = min(dataGraph);
            maxDataGraph = max(dataGraph);
            diffDataGraph = maxDataGraph - minDataGraph;
            %if diffDataGraph == 0 && minDataGraph == 0
            %    ylim(graph, [minDataGraph - 1, maxDataGraph + 1]);
            %else
            %    ylim(graph, [minDataGraph - 0.2*diffDataGraph, maxDataGraph + 0.2*diffDataGraph]);
            %end

            % Row label
            %text(graph(2), mean(xlim(graph(2)))*0.5, mean(ylim(graph(2))), ...
            %    shocksToPlotNames{iShock}, ...
            %    'Interpreter', 'latex', ...
            %    'FontSize', 20);

            % First row
            if iShock == 1
               lg = legend(graph(1), legendStr, 'Location', 'best', 'FontSize', fontSize, 'Orientation', 'horizontal');
               lg.Layout.Tile = 'North';
               %graph(1).Legend.Location = 'North';
            end
            

            % Hide after the first row
            %if iShock > 1
            %    title(graph(1:end), '');
            %    legend(graph(1), 'off');
            %end

        %     % Between lines 1 and 3
        %     if iShock >= 1 && iShock <= 3
        %         set(graph(4:6), 'Color', [0.9 0.9 0.9]);
        %     end

            % Hide after column 2
            %ylabel(graph(4:end), '');

            % Hide after column 2
            %yticks(graph(3:end), []);

            % Turning its axis off
            %axis(graph(2), 'off');

        end

        % Custom font
        set(f.Children.Children, 'FontSize', 12);
        set(f.Children.Children, 'FontWeight', 'bold');

        % Resize title
        for iGraph = 1:length(tl.Children)
            tl.Children(iGraph).Title.FontSize = 16;
        end

        if polRule_simple
           title_polRule_simple = 'Simplified Rule';
        else
           title_polRule_simple = 'Full Rule';
        end

        if stickyPrices
           title_stickyPrices = 'Sticky'; 
        else
           title_stickyPrices = 'Flexible'; 
        end

        set(f, 'Position',  [100, 0, 1400, 150*nShocks]); % resize figure
        exportgraphics(f, ...
            strjoin({pathImages 'IRFs' 'Risky vs Risk-Free' title_stickyPrices ['IRF - Regime ' num2str(regimeId)  ' - All Shocks - ' 'IRF Type ' irf_type ' - PolicyRule ' policyRule ' - ' title_polRule_simple ' - Prices ' title_stickyPrices '.png']}, filesep));
        edu_GraphSetInterpreter(previousInterpreter);
    end
end

%% Plot assets IRFs to different shocks in the same graph - Real Sector (All regimes in the same graph)

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
refMdl = 1; % Pick one model only
mdlGraphVector = mdlVector(refMdl); 
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau     = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 40;
irf_type                                = 'irf'; % 'irf', 'girf'
irf_regime_specific                     = true;
irf_draws                               = 1;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_regime                            = 1:nRegimes;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = nIRFPeriods;
simul_burn                              = 0;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%%%%% PLOT
close all;

varlist = {'polDef', 'y', 'c', 'n', 'bY'};
plotConfig;
modelName   = riseModelNames{refMdl};
defProcess  = riseConfigVector(1).conf_defRProcessType;
policyRule  = riseConfigVector(1).conf_policyRule;

shocksToPlotList        = {};
shocksToPlotNamesList   = {};

shocksToPlotList{1}            = {'epsA', 'epsG', 'epsM'}; % {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'}
shocksToPlotNamesList{1}       = {'$\varepsilon^A$', '$\varepsilon^G$', '$\varepsilon^\mathcal{M}$'};
%shocksToPlotList{2}            = {'epsPolDef'}; % {'epsPolDef', 'epsBeta', 'epsTau'}
%shocksToPlotNamesList{2}       = {'$\varepsilon^{\mathcal{D}}$', '$\varepsilon^\beta$', '$\varepsilon^\tau$'};

plotVars_assets         = {'y', 'c', 'n', 'bY', 'nrPolicy', 'rPolicy', 'Pii'};
plotVarsNames_assets    = {'$Y_t$', '$C_t$', '$N_t$', '${\frac{B}{Y}}_t$', '$i_t$', '$r_t$', '$\Pi_t$'};

nCols = length(plotVars_assets);
nFigs = length(shocksToPlotList);

fontSize = 12;
pColors = linspecer(3);
%customizeLines = false;
customizeLines = struct();

iFigure = 1;
shocksToPlot        = shocksToPlotList{iFigure};
shocksToPlotNames   = shocksToPlotNamesList{iFigure};
nShocks = length(shocksToPlot);

idxAxtiveShocks         = ismember(shocksToPlot, fieldnames(myirfs));
shocksToPlot            = shocksToPlot(idxAxtiveShocks);
shocksToPlotNames       = shocksToPlotNames(idxAxtiveShocks);

% Simplify legend
legendStr_simplified = setUpNames;
legendStr_simplified = replace(legendStr_simplified, ...
                    '$\overline{\mathcal{D}} \approx 0\%$: adjusted by $\mathcal{D}_t$', ...
                    '$\overline{\mathcal{D}} \approx 0\%$');
legendStr_simplified = replace(legendStr_simplified, ...
                    '$\overline{\mathcal{D}} = 5\%$: adjusted by $\mathcal{D}_t$', ...
                    '$\overline{\mathcal{D}} = 5\%$');

f = figure;
set(gca, 'FontSize', 14);
previousInterpreter = edu_GraphSetInterpreter('latex');

% Set layout and reduce empty spaces in the graph
tl = tiledlayout(numel(shocksToPlot), nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';

regimeList = 1:2;
legendStr = [];
for iShock = 1:numel(shocksToPlot)    
    
    legendSubPlot = 1;
    tl_h = nexttile((iShock-1)*nCols + 1);
    for iVar = 2:length(plotVars_assets)
        tl_h = [tl_h, nexttile()];
    end
    
    for regimeId = regimeList
        
        %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);

        if regimeId == regimeList(1)
            customizeLines.custom_LineWidth = {1, 1};
            customizeLines.custom_LineStyle = {'-', '--'};
            customizeLines.custom_LineColor = {'blue', 'red'};
            customizeLines.custom_LineMarker = {'none', 'none'};
        elseif regimeId == regimeList(2)
            customizeLines.custom_LineWidth = {1, 1};
            customizeLines.custom_LineStyle = {':', '-.'};
            customizeLines.custom_LineColor = {'blue', 'red'};
            customizeLines.custom_LineMarker = {'none', 'none'};
        end
        
        if size(endo_ss,1) == 1
            endo_ss_aux = [endo_ss{1,1}]; 
        else
            endo_ss_aux = endo_ss(:,regimeId)';
            endo_ss_aux = [endo_ss_aux{:,:}];
        end
        endo_ss_names = endo_names;
        [irfLevel, irfLevelChg, irfPercChg] = ...
            structIRFsAsContainers(myirfs.(shocksToPlot{iShock}), endo_ss_names, endo_ss_aux, regimeId);

        graph = plotIRFs(tl_h, plotVars_assets, nPlotPeriods, suptitleStr, [], ...
                irfLevel, irfLevelChg, irfPercChg, nCols, plotVarsNames_assets, ...
                0, setFirstPeriodToNaN, customizeLines);

        ylabel(tl_h(1), {shocksToPlotNames{iShock};tl_h(1).YLabel.String}, 'FontSize', 20, 'Interpreter', 'latex');

        % Format subplots
        for iAx = 1:length(graph)
            edu_GraphDrawZeroAxis(graph(iAx));
            box(graph(iAx));
            legend('off');
        end

        % Link axes
        %linkaxes(graph);
        dataGraph = cell2mat(arrayfun(@(x) [x.Children.YData], graph(legendSubPlot:end), 'UniformOutput', false));
        minDataGraph = min(dataGraph);
        maxDataGraph = max(dataGraph);
        diffDataGraph = maxDataGraph - minDataGraph;
        %if diffDataGraph == 0 && minDataGraph == 0
        %    ylim(graph, [minDataGraph - 1, maxDataGraph + 1]);
        %else
        %    ylim(graph, [minDataGraph - 0.2*diffDataGraph, maxDataGraph + 0.2*diffDataGraph]);
        %end

        % Row label
        %text(graph(2), mean(xlim(graph(2)))*0.5, mean(ylim(graph(2))), ...
        %    shocksToPlotNames{iShock}, ...
        %    'Interpreter', 'latex', ...
        %    'FontSize', 20);

        % First row
        if iShock == 1
           legendStr = [legendStr, strcat(['Regime ' num2str(regimeId) ': '], legendStr_simplified)];
        end

         % First row
        if iShock == 1 && regimeId == regimeList(end)
           graphAxes = [graph(1).Children(4:3); graph(1).Children(2:1)];
           lg = legend(graphAxes, legendStr, 'Location', 'best', 'FontSize', fontSize, 'Orientation', 'horizontal');
           lg.Layout.Tile = 'North';
           %graph(1).Legend.Location = 'North';
        end
        

        % Hide after the first row
        %if iShock > 1
        %    title(graph(1:end), '');
        %    legend(graph(1), 'off');
        %end

    %     % Between lines 1 and 3
    %     if iShock >= 1 && iShock <= 3
    %         set(graph(4:6), 'Color', [0.9 0.9 0.9]);
    %     end

        % Hide after column 2
        %ylabel(graph(4:end), '');

        % Hide after column 2
        %yticks(graph(3:end), []);

        % Turning its axis off
        %axis(graph(2), 'off');

    end

    % Custom font
    set(f.Children.Children, 'FontSize', 12);
    set(f.Children.Children, 'FontWeight', 'bold');

    % Resize title
    for iGraph = 1:length(tl.Children)
        tl.Children(iGraph).Title.FontSize = 16;
    end

end


if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

set(f, 'Position',  [100, 0, 1400, 200*nShocks]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'IRFs' modelName 'DefR' defProcess 'RF vs Risky' ['Graph - Real Sector - All Regimes - All Shocks - ' 'IRF Type ' irf_type ' - PolicyRule ' policyRule ' - ' title_polRule_simple ' ' num2str(iFigure) '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Plot assets IRFs to different shocks in the same graph - Real Sector (GIRF)

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlGraphVector = mdlVector;
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau     = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 40;
irf_type                                = 'girf'; % 'irf', 'girf'
irf_regime_specific                     = false;
irf_draws                               = 1000;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_regime                            = 1:4;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = nIRFPeriods;
simul_burn                              = 100;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false

% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%%%%% PLOT
close all;

refMdl = 1;
varlist = {'polDef', 'y', 'c', 'n', 'bY'};
plotConfig;
modelName   = riseModelNames{refMdl};
defProcess  = riseConfigVector(1).conf_defRProcessType;
policyRule  = riseConfigVector(1).conf_policyRule;

shocksToPlotList        = {};
shocksToPlotNamesList   = {};

shocksToPlotList{1}            = {'epsA', 'epsG', 'epsM'}; % {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'}
shocksToPlotNamesList{1}       = {'$\varepsilon^A$', '$\varepsilon^G$', '$\varepsilon^\mathcal{M}$'};
%shocksToPlotList{2}            = {'epsPolDef'}; % {'epsPolDef', 'epsBeta', 'epsTau'}
%shocksToPlotNamesList{2}       = {'$\varepsilon^{\mathcal{D}}$', '$\varepsilon^\beta$', '$\varepsilon^\tau$'};

nPlotPeriods            = 1:nIRFPeriods;
plotVars_assets         = {'y', 'c', 'n', 'bY', 'nrPolicy', 'rPolicy', 'Pii'};
plotVarsNames_assets    = {'$Y_t$', '$C_t$', '$N_t$', '${\frac{B}{Y}}_t$', '$i_t$', '$r_t$', '$\Pi_t$'};

nCols = length(plotVars_assets);
nFigs = length(shocksToPlotList);

fontSize = 12;
pColors = linspecer(3); 
%customizeLines = false;
%customizeLines = struct();
customizeLines.custom_LineWidth = {1, 2, 1, 1, 1, 1, 1};
customizeLines.custom_LineStyle = {'-', ':', '--', '-.', '-', ':', '--'};
customizeLines.custom_LineColor = {pColors(1,:), pColors(2,:), pColors(3,:), 'blue', 'blue', 'blue', 'blue'};
customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};

if irf_regime_specific == false
    nReg = 1;
else
    nReg = 4;
end
for regimeId = 1:nReg
    for iFigure = 1:nFigs

        shocksToPlot        = shocksToPlotList{iFigure};
        shocksToPlotNames   = shocksToPlotNamesList{iFigure};
        nShocks = length(shocksToPlot);

        idxAxtiveShocks         = ismember(shocksToPlot, fieldnames(myirfs));
        shocksToPlot            = shocksToPlot(idxAxtiveShocks);
        shocksToPlotNames       = shocksToPlotNames(idxAxtiveShocks);

        f = figure;
        set(gca, 'FontSize', 14);
        previousInterpreter = edu_GraphSetInterpreter('latex');

        % Set layout and reduce empty spaces in the graph
        tl = tiledlayout(numel(shocksToPlot), nCols);
        tl.TileSpacing  = 'compact';
        tl.Padding      = 'compact';
        
        for iShock = 1:numel(shocksToPlot)

            %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);

            if size(endo_ss,1) == 1
                endo_ss_aux = [endo_ss{1,1}]; 
            else
                endo_ss_aux = endo_ss(:,regimeId)';
                endo_ss_aux = [endo_ss_aux{:,:}];
            end
            endo_ss_names = endo_names;
            [irfLevel, irfLevelChg, irfPercChg] = ...
                structIRFsAsContainers(myirfs.(shocksToPlot{iShock}), endo_ss_names, endo_ss_aux, regimeId);

            legendSubPlot = 1;
            tl_h = nexttile();
            for iVar = 2:length(plotVars_assets)
                tl_h = [tl_h, nexttile()];
            end

            graph = plotIRFs(tl_h, plotVars_assets, nPlotPeriods, suptitleStr, [], ...
                    irfLevel, irfLevelChg, irfPercChg, nCols, plotVarsNames_assets, ...
                    0, setFirstPeriodToNaN, customizeLines);
            
            ylabel(tl_h(1), {shocksToPlotNames{iShock};tl_h(1).YLabel.String}, 'FontSize', 20, 'Interpreter', 'latex');
                
            % Format subplots
            for iAx = 1:length(graph)
                edu_GraphDrawZeroAxis(graph(iAx));
                box(graph(iAx));
                legend('off');
            end

            % Link axes
            %linkaxes(graph);
            dataGraph = cell2mat(arrayfun(@(x) [x.Children.YData], graph(legendSubPlot:end), 'UniformOutput', false));
            minDataGraph = min(dataGraph);
            maxDataGraph = max(dataGraph);
            diffDataGraph = maxDataGraph - minDataGraph;
            %if diffDataGraph == 0 && minDataGraph == 0
            %    ylim(graph, [minDataGraph - 1, maxDataGraph + 1]);
            %else
            %    ylim(graph, [minDataGraph - 0.2*diffDataGraph, maxDataGraph + 0.2*diffDataGraph]);
            %end

            % Row label
            %text(graph(2), mean(xlim(graph(2)))*0.5, mean(ylim(graph(2))), ...
            %    shocksToPlotNames{iShock}, ...
            %    'Interpreter', 'latex', ...
            %    'FontSize', 20);

            % First row
            if iShock == 1
               lg = legend(graph(1), legendStr, 'Location', 'best', 'FontSize', fontSize, 'Orientation', 'horizontal');
               lg.Layout.Tile = 'North';
               %graph(1).Legend.Location = 'North';
            end
            

            % Hide after the first row
            %if iShock > 1
            %    title(graph(1:end), '');
            %    legend(graph(1), 'off');
            %end

        %     % Between lines 1 and 3
        %     if iShock >= 1 && iShock <= 3
        %         set(graph(4:6), 'Color', [0.9 0.9 0.9]);
        %     end

            % Hide after column 2
            %ylabel(graph(4:end), '');

            % Hide after column 2
            %yticks(graph(3:end), []);

            % Turning its axis off
            %axis(graph(2), 'off');

        end

        % Custom font
        set(f.Children.Children, 'FontSize', 12);
        set(f.Children.Children, 'FontWeight', 'bold');

        % Resize title
        for iGraph = 1:length(tl.Children)
            tl.Children(iGraph).Title.FontSize = 16;
        end

        if polRule_simple
           title_polRule_simple = 'Simplified Rule';
        else
           title_polRule_simple = 'Full Rule';
        end

        if irf_regime_specific == false
            regimeStr = 'RSS';
        else
            regimeStr = num2str(regimeId);
        end
        

        set(f, 'Position',  [100, 0, 1400, 220*nShocks]); % resize figure
        exportgraphics(f, ...
            strjoin({pathImages 'IRFs' modelName 'DefR' defProcess 'RF vs Risky' ['Graph - Real Sector - Regime ' regimeStr  ' - All Shocks - ' 'IRF Type ' irf_type ' - PolicyRule ' policyRule ' - ' title_polRule_simple ' ' num2str(iFigure) '.png']}, filesep));
        edu_GraphSetInterpreter(previousInterpreter);
    end
end

%% Plot assets IRFs to different shocks in the same graph - Inflation - All Policy Rules (Regime-Specific)

%%%% SELECT THE REGIMES YOU WANT TO PLOT
regIRFs = 1:2; % 1:2 or 3:4

%%%% UNCOMMENT THE VARIABLE YOU WANT TO PLOT
%plotVars_assets         = {'nrPolicy'};
%plotVarsNames_assets    = {'$i_t$'};
plotVars_assets         = {'Pii'};
plotVarsNames_assets    = {'$\Pi_t$'};

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlGraphVector = mdlVector;
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau     = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 20;
irf_type                                = 'irf'; % 'irf', 'girf'
irf_regime_specific                     = true;
irf_draws                               = 1;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_regime                            = 1:4;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = nIRFPeriods;
simul_burn                              = 0;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false


% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%%%%% PLOT
%close all;

refMdl = 1;
varlist = {'polDef', 'y', 'c', 'n', 'bY'};
plotConfig;
modelName   = riseModelNames{refMdl};
defProcess  = riseConfigVector(1).conf_defRProcessType;
policyRule  = riseConfigVector(1).conf_policyRule;
stickyPrices = conf_stickyPrices;

shocksToPlotList        = {};
shocksToPlotNamesList   = {};

shocksToPlotList{1}            = {'epsA', 'epsG', 'epsM'}; % {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'}
shocksToPlotNamesList{1}       = {'$\varepsilon^A$', '$\varepsilon^G$', '$\varepsilon^\mathcal{M}$'};
%shocksToPlotList{2}            = {'epsPolDef'}; % {'epsPolDef', 'epsBeta', 'epsTau'}
%shocksToPlotNamesList{2}       = {'$\varepsilon^{\mathcal{D}}$', '$\varepsilon^\beta$', '$\varepsilon^\tau$'};

nReg = length(regIRFs);
nVars = length(plotVars_assets);
nCols = nReg * nVars;
nFigs = length(shocksToPlotList);

nPlotPeriods = 1:nIRFPeriods;

fontSize = 12;
pColors = linspecer(3); 
%customizeLines = false;
%customizeLines = struct();
N = 4;
C = linspecer(N);
customizeLines.custom_LineWidth = {2, 2, 2, 2, 2, 2, 2};
customizeLines.custom_LineStyle = {'-', '-', '-', '-', '-', ':', '-.'};
customizeLines.custom_LineColor = {C(1,:), C(2,:), C(3,:), C(4,:), 'blue', 'blue', 'blue'};
customizeLines.custom_LineMarker = {'none', 'o', '+', 's', 'none', 'none', 'none'};

f = figure;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'FontSize', 14);
previousInterpreter = edu_GraphSetInterpreter('latex');

iFigure = 1;

shocksToPlot        = shocksToPlotList{iFigure};
shocksToPlotNames   = shocksToPlotNamesList{iFigure};
nShocks = length(shocksToPlot);

idxAxtiveShocks         = ismember(shocksToPlot, fieldnames(myirfs));
shocksToPlot            = shocksToPlot(idxAxtiveShocks);
shocksToPlotNames       = shocksToPlotNames(idxAxtiveShocks);

% Set layout and reduce empty spaces in the graph
tl = tiledlayout(numel(shocksToPlot), nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';
        
for regimeId = regIRFs

        %title(tl, 'Regime 1 \quad\quad\quad\quad\quad\quad\quad\quad\quad Regime 2', 'FontSize', 20, 'Interpreter', 'latex');
    
        for iShock = 1:numel(shocksToPlot)

            %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);

            if size(endo_ss,1) == 1
                endo_ss_aux = [endo_ss{1,1}]; 
            else
                endo_ss_aux = endo_ss(:,regimeId)';
                endo_ss_aux = [endo_ss_aux{:,:}];
            end
            endo_ss_names = endo_names;
            [irfLevel, irfLevelChg, irfPercChg] = ...
                structIRFsAsContainers(myirfs.(shocksToPlot{iShock}), endo_ss_names, endo_ss_aux, regimeId);

            legendSubPlot = 1;
            tl_h = [nexttile((iShock-1)*nCols + (regimeId-regIRFs(1))*nVars + 1)];
            for iVar = 2:length(plotVars_assets)
                tl_h = [tl_h, nexttile()];
            end

            graph = plotIRFs(tl_h, plotVars_assets, nPlotPeriods, suptitleStr, [], ...
                    irfLevel, irfLevelChg, irfPercChg, 1, plotVarsNames_assets, ...
                    0, setFirstPeriodToNaN, customizeLines);
            
            % Replace nrPolicy by the difference w.r.t. specific rule
%             ruleRef = 2;
%             graph(2).Children(1).YData = graph(2).Children(1).YData - graph(2).Children(ruleRef).YData;
%             graph(2).Children(3).YData = graph(2).Children(3).YData - graph(2).Children(ruleRef).YData;
%             graph(2).Children(2).YData = graph(2).Children(2).YData - graph(2).Children(ruleRef).YData;
%             % Adjust y axis
%             ax = graph(2);
%             graphGapPct = 0.1;
%             yData = horzcat(ax.Children.YData); % concatenate all series
%             if min([yData, 0]) == max([yData, 0])
%                 ax.YLim = [-1 1];
%             elseif max(abs(yData)) < tol
%                 ax.YLim = [-1 1];
%             elseif sum(isfinite(max(abs(yData)))) == 0
%                 ax.YLim = [-1 1];
%             else
%                 ax.YLim = [min([yData, 0]) max([yData, 0])] * (1+graphGapPct);
%             end
%             %%%%%%%%%%%%%%%%
            
            % Row shock names 
            if regimeId == 1 || regimeId == 3
                ylabel(tl_h(1), {shocksToPlotNames{iShock};tl_h(1).YLabel.String}, 'FontSize', 20, 'Interpreter', 'latex');
            end
            
            % Column var name
            if iShock == 1
                colTitles = {['Regime ' num2str(regimeId)], plotVarsNames_assets{1}}';
                title(tl_h, colTitles, 'FontSize', 20, 'Interpreter', 'latex');
            else
                title(tl_h, '');
            end
                
            % Format subplots
            for iAx = 1:length(graph)
                edu_GraphDrawZeroAxis(graph(iAx));
                box(graph(iAx));
                legend('off');
                graph(iAx).XGrid = 'on';
            end

            % Link axes
            %linkaxes(graph);
            dataGraph = cell2mat(arrayfun(@(x) [x.Children.YData], graph(legendSubPlot:end), 'UniformOutput', false));
            minDataGraph = min(dataGraph);
            maxDataGraph = max(dataGraph);
            diffDataGraph = maxDataGraph - minDataGraph;
            %if diffDataGraph == 0 && minDataGraph == 0
            %    ylim(graph, [minDataGraph - 1, maxDataGraph + 1]);
            %else
            %    ylim(graph, [minDataGraph - 0.2*diffDataGraph, maxDataGraph + 0.2*diffDataGraph]);
            %end

            % First row
            if iShock == 1 && regimeId == regIRFs(1)
               lg = legend(graph(1), legendStr, 'Location', 'best', 'FontSize', fontSize, 'Orientation', 'horizontal');
               lg.Layout.Tile = 'North';
               %graph(1).Legend.Location = 'North';
            end

        end   

end

% Custom font
set(f.Children.Children, 'FontSize', 12);
set(f.Children.Children, 'FontWeight', 'bold');

% Resize title
for iGraph = 1:length(tl.Children)
    tl.Children(iGraph).Title.FontSize = 16;
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

if stickyPrices
   title_stickyPrices = 'Sticky'; 
else
   title_stickyPrices = 'Flexible'; 
end

set(f, 'Position',  [100, 0, 1000, 200*nShocks]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'IRFs' 'Rules Comparison' title_stickyPrices ['Graph - Inflation - Regimes ' num2str(regIRFs(1)) num2str(regIRFs(2)) ' - All Shocks - ' 'IRF Type ' irf_type ' - PolicyRule All - ' title_polRule_simple ' - Prices ' title_stickyPrices ' - Variable ' plotVars_assets{1} '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% Plot assets IRFs to different shocks in the same graph - Inflation - All Policy Rules (GIRF)

polRule_simple = true;
selectedShocksTurnOff = true;

%%%%%%%%%%%%%% SET-UP MODELS
mdlGraphVector = mdlVector;
for iMdl = 1:length(mdlGraphVector)
    paramsStructTemp = struct();
    if polRule_simple
        %paramsStructTemp.piiBar  = 0.011 ;
        paramsStructTemp.phi_dY  = 0 ;
        paramsStructTemp.phi_Y   = 0 ;
        paramsStructTemp.phiNr   = 0 ;
        paramsStructTemp.phi_Exp = 0 ;
    end

    % Turn off shocks
    if selectedShocksTurnOff
        paramsStructTemp.sigmaTau     = 0 ;
        paramsStructTemp.sigmaBeta    = 0 ;
        %paramsStructTemp.sigmaPolDef = 0 ;
    end

    mdlGraphVector(iMdl)      = set(mdlGraphVector(iMdl),'parameters', paramsStructTemp);
    mdlGraphVector(iMdl)      = solve(mdlGraphVector(iMdl));
end
%%%%%%%%%%%%%% COMPUTE IRFs

% Initiate randomizer
rng(1900)

useParallel = true;

shockSize = 1;
%M=[set(mdlVector,  'irf_shock_sign', -shockSize),...
%    set(mdlVector, 'irf_shock_sign',  shockSize)];
M = set(mdlGraphVector,  'irf_shock_sign', shockSize);

nIRFPeriods                             = 20;
irf_type                                = 'girf'; % 'irf', 'girf'
irf_regime_specific                     = true;
irf_draws                               = 10000;
irf_girf_regime_uncertainty             = true; % If "false" then impose the sequence of regimes from the first simulation to the second
simul_regime                            = 1:4;
simul_order                             = 1;
simul_pruned                            = true;
simul_honor_constraints                 = true; % true or false
simul_honor_constraints_through_switch  = true; % true or false
simul_anticipate_zero                   = true;
simul_periods                           = nIRFPeriods;
simul_burn                              = 1000;
simul_frwrd_back_shoot                  = false; % true or false
simul_shock_uncertainty                 = true; % true or false


% IRF is calculated as absolute deviation from the steady state (or the imposed solution point)
tic;
if useParallel
    if isempty(gcp('nocreate'))
        parpool; % Start parallelization
    end
end
computeIRFs;
toc;

%%%%%%%%%%%%%%%%% PLOT
%close all;

refMdl = 1;
varlist = {'polDef', 'y', 'c', 'n', 'bY'};
plotConfig;
modelName   = riseModelNames{refMdl};
defProcess  = riseConfigVector(1).conf_defRProcessType;
policyRule  = riseConfigVector(1).conf_policyRule;

shocksToPlotList        = {};
shocksToPlotNamesList   = {};

shocksToPlotList{1}            = {'epsA', 'epsG', 'epsM'}; % {'epsA', 'epsR', 'epsD', 'epsG', 'epsM'}
shocksToPlotNamesList{1}       = {'$\varepsilon^A$', '$\varepsilon^G$', '$\varepsilon^\mathcal{M}$'};
%shocksToPlotList{2}            = {'epsPolDef'}; % {'epsPolDef', 'epsBeta', 'epsTau'}
%shocksToPlotNamesList{2}       = {'$\varepsilon^{\mathcal{D}}$', '$\varepsilon^\beta$', '$\varepsilon^\tau$'};

plotVars_assets         = {'Pii'};
plotVarsNames_assets    = {'$\Pi_t$'};

regIRFs = 1:2; % 1:2 or 3:4
nReg = length(regIRFs);
nVars = length(plotVars_assets);
nCols = nReg * nVars;
nFigs = length(shocksToPlotList);

nPlotPeriods = 1:nIRFPeriods;

fontSize = 12;
pColors = linspecer(3); 
%customizeLines = false;
%customizeLines = struct();
customizeLines.custom_LineWidth = {1, 2, 1, 1, 1, 1, 1};
customizeLines.custom_LineStyle = {'-', ':', '-.', '--', '-', ':', '-.'};
customizeLines.custom_LineColor = {'blue', 'red', 'black', 'blue', 'blue', 'blue', 'blue'};
customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};

f = figure;
set(gca, 'FontSize', 14);
previousInterpreter = edu_GraphSetInterpreter('latex');

iFigure = 1;

shocksToPlot        = shocksToPlotList{iFigure};
shocksToPlotNames   = shocksToPlotNamesList{iFigure};
nShocks = length(shocksToPlot);

idxAxtiveShocks         = ismember(shocksToPlot, fieldnames(myirfs));
shocksToPlot            = shocksToPlot(idxAxtiveShocks);
shocksToPlotNames       = shocksToPlotNames(idxAxtiveShocks);

% Set layout and reduce empty spaces in the graph
tl = tiledlayout(numel(shocksToPlot), nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';
        
for regimeId = regIRFs

        %title(tl, 'Regime 1 \quad\quad\quad\quad\quad\quad\quad\quad\quad Regime 2', 'FontSize', 20, 'Interpreter', 'latex');
    
        for iShock = 1:numel(shocksToPlot)

            %graph = figure('name',['Impulse responses to a ' shock_names{iShock} 'shock']);

            if size(endo_ss,1) == 1
                endo_ss_aux = [endo_ss{1,1}]; 
            else
                endo_ss_aux = endo_ss(:,regimeId)';
                endo_ss_aux = [endo_ss_aux{:,:}];
            end
            endo_ss_names = endo_names;
            [irfLevel, irfLevelChg, irfPercChg] = ...
                structIRFsAsContainers(myirfs.(shocksToPlot{iShock}), endo_ss_names, endo_ss_aux, regimeId);

            legendSubPlot = 1;
            tl_h = [nexttile((iShock-1)*nCols + (regimeId-regIRFs(1))*nVars + 1)];
            for iVar = 2:length(plotVars_assets)
                tl_h = [tl_h, nexttile()];
            end

            graph = plotIRFs(tl_h, plotVars_assets, nPlotPeriods, suptitleStr, [], ...
                    irfLevel, irfLevelChg, irfPercChg, 1, plotVarsNames_assets, ...
                    0, setFirstPeriodToNaN, customizeLines);
            
            % Replace nrPolicy by the difference w.r.t. specific rule
%             ruleRef = 2;
%             graph(2).Children(1).YData = graph(2).Children(1).YData - graph(2).Children(ruleRef).YData;
%             graph(2).Children(3).YData = graph(2).Children(3).YData - graph(2).Children(ruleRef).YData;
%             graph(2).Children(2).YData = graph(2).Children(2).YData - graph(2).Children(ruleRef).YData;
%             % Adjust y axis
%             ax = graph(2);
%             graphGapPct = 0.1;
%             yData = horzcat(ax.Children.YData); % concatenate all series
%             if min([yData, 0]) == max([yData, 0])
%                 ax.YLim = [-1 1];
%             elseif max(abs(yData)) < tol
%                 ax.YLim = [-1 1];
%             elseif sum(isfinite(max(abs(yData)))) == 0
%                 ax.YLim = [-1 1];
%             else
%                 ax.YLim = [min([yData, 0]) max([yData, 0])] * (1+graphGapPct);
%             end
%             %%%%%%%%%%%%%%%%
            
            % Row shock names 
            if regimeId == 1
                ylabel(tl_h(1), {shocksToPlotNames{iShock};tl_h(1).YLabel.String}, 'FontSize', 20, 'Interpreter', 'latex');
            end
            
            % Column var name
            if iShock == 1
                colTitles = {['Regime ' num2str(regimeId)], plotVarsNames_assets{1}}';
                title(tl_h, colTitles, 'FontSize', 20, 'Interpreter', 'latex');
            else
                title(tl_h, '');
            end
                
            % Format subplots
            for iAx = 1:length(graph)
                edu_GraphDrawZeroAxis(graph(iAx));
                box(graph(iAx));
                legend('off');
            end

            % Link axes
            %linkaxes(graph);
            dataGraph = cell2mat(arrayfun(@(x) [x.Children.YData], graph(legendSubPlot:end), 'UniformOutput', false));
            minDataGraph = min(dataGraph);
            maxDataGraph = max(dataGraph);
            diffDataGraph = maxDataGraph - minDataGraph;
            %if diffDataGraph == 0 && minDataGraph == 0
            %    ylim(graph, [minDataGraph - 1, maxDataGraph + 1]);
            %else
            %    ylim(graph, [minDataGraph - 0.2*diffDataGraph, maxDataGraph + 0.2*diffDataGraph]);
            %end

            % Row label
            %text(graph(2), mean(xlim(graph(2)))*0.5, mean(ylim(graph(2))), ...
            %    shocksToPlotNames{iShock}, ...
            %    'Interpreter', 'latex', ...
            %    'FontSize', 20);

            % First row
            if iShock == 1 && regimeId == 1
               lg = legend(graph(1), legendStr, 'Location', 'best', 'FontSize', fontSize, 'Orientation', 'horizontal');
               lg.Layout.Tile = 'North';
               %graph(1).Legend.Location = 'North';
            end
            

            % Hide after the first row
            %if iShock > 1
            %    title(graph(1:end), '');
            %    legend(graph(1), 'off');
            %end

        %     % Between lines 1 and 3
        %     if iShock >= 1 && iShock <= 3
        %         set(graph(4:6), 'Color', [0.9 0.9 0.9]);
        %     end

            % Hide after column 2
            %ylabel(graph(4:end), '');

            % Hide after column 2
            %yticks(graph(3:end), []);

            % Turning its axis off
            %axis(graph(2), 'off');

        end   

end

% Custom font
set(f.Children.Children, 'FontSize', 12);
set(f.Children.Children, 'FontWeight', 'bold');

% Resize title
for iGraph = 1:length(tl.Children)
    tl.Children(iGraph).Title.FontSize = 16;
end

if polRule_simple
   title_polRule_simple = 'Simplified Rule';
else
   title_polRule_simple = 'Full Rule';
end

set(f, 'Position',  [100, 0, 800, 200*nShocks]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'IRFs' 'Rules Comparison' title_stickyPrices ['Graph - Inflation - Regimes ' num2str(regIRFs(1)) num2str(regIRFs(2)) ' - All Shocks - ' 'IRF Type ' irf_type ' - PolicyRule All - ' title_polRule_simple ' - Prices ' title_stickyPrices ' - Variable ' plotVars_assets{1} '.png']}, filesep));
edu_GraphSetInterpreter(previousInterpreter);

%% GIRF
%%%%% NÃO TÁ FUNCIONANDO
girf_draws = 1000;
simul_burn = 200;
simul_periods = 20;
diff_y = NaN(girf_draws,1);
for irf_iter = 1:girf_draws
    %generate baseline shocks and simulate
    sims_baseline = simulate(mdlVector, 'simul_periods', simul_periods, 'simul_burn', simul_burn);
    %add deterministic impulse and simulate
    sims_baseline.epsA(1,:) = sims_baseline.epsA(1,:) + 1;
    sims_shock = simulate(mdlVector,  'simul_periods', simul_periods, 'simul_burn', simul_burn, 'simul_historical_data', sims_baseline, 'cond_exo_vars',{'epsA','epsG','epsM'});
    % obtain the difference
    diff_y(irf_iter) = sims_shock.y - sims_baseline.y;
end
