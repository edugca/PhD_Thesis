%% Load parameters / Simulate fiscal limit

% Should recalculate the fiscal limits?
createFiscalLimitPolicyFunctions    = false;
saveFiscalLimitPolicyFunctions      = false;
loadFiscalLimitPolicyFunctions      = true;

% Define language
useLanguage = 'MatLab'; % 'C++' or 'MatLab' or 'GPU'
if strcmp(useLanguage, 'GPU')
    numVarType = 'single';
else
    numVarType = 'single';
end

policyRule          = riseConfig.conf_policyRule;
hasGovernment       = riseConfig.conf_hasGovernment;
govExpenses         = riseConfig.conf_govExpenses;
debtLevel           = riseConfig.conf_debtLevel;
defPolTarget        = riseConfig.conf_defPolTarget;
govAccumDebt        = riseConfig.conf_govAccumDebt;
govBondsAreRiskFree = riseConfig.conf_govBondsAreRiskFree;
defRProcessType     = riseConfig.conf_defRProcessType;
prodProcessType     = riseConfig.conf_prodProcessType;
occBinConstraint    = riseConfig.conf_occBinConstraint;
calibrationName     = riseConfig.conf_calibrationName;
shocksTurnedOff     = riseConfig.conf_shocksTurnedOff;
simplifyPolicyRule  = riseConfig.conf_simplifyPolicyRule;
approximateDefaultProb = riseConfig.conf_approximateDefaultProb;
stickyPrices        = riseConfig.conf_stickyPrices;
debtLevelGuess      = riseConfig.conf_debtLevelGuess;

nPeriods        = 1;          % t = 1 is Steady State
nFiscalLimit    = 200;        % Number of periods at each simulation of the fiscal limit
nSimulations    = 50000;      % Number of simulations of the fiscal limit
nSimIRFs        = 1;          % Number of IRFs to be simulated
nBurnIn         = 200;        % Number of periods of shocks to be excluded
tol             = 1e-6;       % Tolerance for convergence


simStruct = struct();
load('shocks.mat');
load('fiscalLimitPolicyFunction.mat');

if paperNumber == 2
    paper2_params_loadStandard;        % Parameter Values
elseif paperNumber == 3
    paper3_params_loadStandard;        % Parameter Values
end
params_loadValues;

%approximateFiscalLimits;

ssValues_loadStandard;  % Steady-state values
ssValues_loadValues;

if createFiscalLimitPolicyFunctions

    % Simulate the fiscal limits
    simulate_fiscalLimit;
    
    % Load parameters and steady state values
    load('shocks.mat');
    load('fiscalLimitPolicyFunction.mat');
    run('params_loadStandard.m');        % Parameter Values
    run('params_loadValues.m');
    run('ssValues_loadStandard.m');  % Steady-state values
    run('ssValues_loadValues.m');
    
end