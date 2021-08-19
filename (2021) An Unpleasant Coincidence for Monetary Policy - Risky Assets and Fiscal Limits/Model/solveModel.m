%% SOLVE MODEL

%%%%%%%%%%% Set up

% Set-Ups
if paperNumber == 2
    paper2_setUps_polRules;
elseif paperNumber == 3
    paper3_setUps_polRules;
end
nMdl        = length(riseConfigVector);
setUpNames  = riseSetUpNames;

for iSetUp = 1:length(riseConfigVector)
   riseConfigVector(iSetUp).conf_calibrationName = conf_calibrationName; 
end

clear('riseMdlVector');
for iMdl=1:nMdl
    
    % Set configuration
    riseConfig = riseConfigVector(iMdl);
    
    % Load parameters / Simulate fiscal limit
    loadParamsAndSimulFiscalLimits;

%     % Approximate fiscal limits by a polynomial function
%     if riseConfigVector(iMdl).conf_hasGovernment
%         approximateFiscalLimits;
%     end
    
    % RISE file path
    if paperNumber == 2
        modelPath = 'paper2_model';
    elseif paperNumber == 3
        modelPath = 'paper3_model';
    end
    
    % Load model
     if riseConfig.conf_estimateModel == false
        mdlVector(iMdl) = rise(modelPath, ...
            'rise_flags', riseConfigVector(iMdl), ...
            ...%'solver', 'mfi', ...            % solvers={'mnk','fwz','mn','mfi'};
            'solve_derivatives_type', solve_derivatives_type, ... % symbolic, automatic, numerical
            'max_deriv_order', solveOrder ...
            ); % symbolic, automatic, numerical
     else    
        mdlVector(iMdl) = rise(modelPath, ...
            'rise_flags', riseConfigVector(iMdl), ...
            ...%'solver', 'mfi', ...            % solvers={'mnk','fwz','mn','mfi'};
            'solve_derivatives_type', solve_derivatives_type, ... % symbolic, automatic, numerical
            'max_deriv_order', solveOrder, ...
            'data', myData ...
            );
    end

    % Add customized options
    mdlVector(iMdl).user_data   = riseConfigVector(iMdl);
    
    % Set alternative calibration
    %paramsStruct = simStruct.params;
    setInitialPriorValuesToParams = true;
    [paramsStruct, priorsStruct] = createParameters(...
        mdlVector(iMdl), riseConfigVector(iMdl), ...
        calibrationName, setInitialPriorValuesToParams, simStruct);
    
    if strcmp(riseConfig.conf_occBinConstraint, 'None')
        % [], or this is the ID of the reference regime);
        solve_occbin = [];
    elseif strcmp(riseConfig.conf_occBinConstraint, 'FiscalLimit')
        
        % assign each restriction to a markov chain
        restr_map = struct('r1taxLim', 1, 'r2fisLim', 2); % struct('taxLim', 1, 'fisLim', 2);
        
%         paramsStruct.gammaTau_r1taxLim_1   = 0.03;
%         paramsStruct.gammaTau_r1taxLim_2   = 0.03;
%         paramsStruct.phi_r1taxLim_1       = 1.5;
%         paramsStruct.phi_r1taxLim_2       = 1.5;
        
        % [], or this is the ID of the reference regime);
        %solve_occbin = {1, restr_map, mdlVector(iMdl).markov_chains.regimes};
        %solve_occbin = 1;
        
        
        solve_occbin = [];
        paramsStruct = rmfield(paramsStruct, 'r1fisLim_tp_1_2');
        paramsStruct = rmfield(paramsStruct, 'r1fisLim_tp_2_1');
        paramsStruct = rmfield(paramsStruct, 'r2taxLim_tp_1_2');
        paramsStruct = rmfield(paramsStruct, 'r2taxLim_tp_2_1');
    end
    
    % Specify parameters of the model
    mdlVector(iMdl) = set(mdlVector(iMdl),'parameters', paramsStruct);
    
    %%% JUST IN CASE ONE NEEDS BETTER UNDERSTANDING OF RISE
    helpRISE;
    %%%%%%%%%%%%%%
    mdlVector(iMdl) = solve(mdlVector(iMdl), ...
        'debug', debugMode, ...
        'steady_state_file', 'paper2_steadyStateFile', ...
        'solve_occbin', solve_occbin, ... % [], or this is the ID of the reference regime);
        'solve_order', solveOrder, ...           % may not work under endogenous regime switching
        ...%'steady_state_fixed',   true, ...  % some variables are fixed
        'steady_state_imposed', steady_state_imposed, ...   % if imposed, then does not check whether it is a steady state
        'steady_state_unique',  steady_state_unique);      % the steady state is unique: solves around the ergodic mean and does not run the steady state file
end

% Auxiliary variables
if strcmp(riseConfig.conf_occBinConstraint, 'FiscalLimit')
    regimeNames     = { 'Tax rate below max / Debt below fiscal limit', ...
                        'Tax rate at max / Debt below fiscal limit', ...
                        'Tax rate below max / Debt crossed fiscal limit', ...
                        'Tax rate at max / Debt crossed fiscal limit'};
elseif strcmp(riseConfig.conf_occBinConstraint, 'None') && strcmp(riseConfig.conf_defRProcessType, 'MS_Exogenous')
    regimeNames     = { 'High growth', ...
                        'Low growth'};
end

endo_names      = arrayfun(@(c) c.endogenous.name, mdlVector, 'UniformOutput', false);
endo_texnames   = arrayfun(@(c) c.endogenous.tex_name, mdlVector, 'UniformOutput', false);
shock_names     = arrayfun(@(c) c.exogenous.name, mdlVector, 'UniformOutput', false);
shock_texnames  = arrayfun(@(c) c.exogenous.tex_name, mdlVector, 'UniformOutput', false);
shock_list      = shock_names;

% Steady states
nRegimes = [];
for iMdl = 1:nMdl
    nRegimes(iMdl) = mdlVector(iMdl).markov_chains.regimes_number;
end

% Steady states
endo_ss     = cell(nMdl, max(nRegimes));
for iMdl = 1:nMdl
    for iReg = 1:nRegimes(iMdl)
        endo_ss{iMdl, iReg} = mdlVector(iMdl).solution.ss{iReg};
    end
end
