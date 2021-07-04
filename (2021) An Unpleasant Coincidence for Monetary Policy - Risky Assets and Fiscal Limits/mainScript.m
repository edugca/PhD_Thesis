%%%%%%%%% Amaral(2021): An Unpleasant Coincidence for Monetary Policy: Risky Assets and Fiscal Limits

%% RUN THIS SECTION BEFORE RUNNING ANY OF THE OTHERS
addpath(genpath('Auxiliary functions'));

addpath('Laffer Curve');
addpath('Calibration');
addpath('Fiscal Limits');
addpath(genpath('Model'));
addpath('Saved');

%% Plot Figure 1
lafferCurve;

%% Plot Figure 2

% Open routine setUps_polRules
% Uncomment the code of any model in the section "GRAPH: Risky 4 Rules: Actual debt level", leaving the other models in that section commented
% Comment the code in other sections
% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'symbolic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate the model"
% Run section "Simulation: PLOT tau_t vs. B_t"
% If the graph plotted does not contain a period in which the peak of the Laffer curve binds, change the start and end periods in rngPeriods = [start end];

%% Calibration
% TFP Loss                              = TFP Loss.xlsx
% Steady-state debt level               = Debt to GDP.xlsx
% Comparison between gross and net debt = Debt and Taxes.xlsx

%% Plot Figure 3: Goodness of Fit

% Open routine setUps_polRules
% Uncomment the code of any model in the section "GRAPH: Risky 4 Rules: Actual debt level", leaving the other models in that section commented
% Comment the code in other sections
% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% Run section "Goodness of fit of the logistic approximation"

%% Plot Figures 4 and 5 and possible calibration of GammaTau
% dataset = GFSMAB_06-13-2020 16-01-20-13_timeSeries.csv
loadGovernmentExpenses;

%% Plot Figures 6 to 10 and Table 1 (LaTex)
% dataset = GFSMAB_06-13-2020 16-01-20-13_timeSeries.csv
distFiscalLimits;

%% To reestimate the fiscal limits

% To reestimate the fiscal limits, follow the instructions in the top
% section of the distFiscalLimits routine

%% Plot Figures 11 to 14 and 39 to 42

% Open routine setUps_polRules
% Uncomment the code in the section "GRAPH: ParamStability: Risky All Rules: Actual debt level"
% Uncomment the code in the section "GRAPH: ParamStability: Risky All Rules: High debt level"
% Comment the code in other sections
% Open mainScript_model
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Test stability"

% It shall take a long time!

%% Plot Figures 15, 16, and 43 to 46

% Open routine setUps_polRules
% Uncomment the code in the section "GRAPH: Risk-free vs. Risky (Rule: target rPolicy)"
% Comment the code in other sections
% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'symbolic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Plot assets IRFs to different shocks in the same graph - Real Sector (Regime-Specific)"

%% Plot Figures 17, 18, and 45 to 46

% Open routine setUps_polRules
% Uncomment the code in the section "GRAPH: Risky 4 Rules: 1% debt level"
% Comment the code in other sections
% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'symbolic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% In the section "Plot assets IRFs to different shocks in the same graph - Inflation - All Policy Rules (Regime-Specific)",
% modify the following parameters right at the beginning of the section:
%
% 1) Select the regimes to plot
% regIRFs = 3:4; % 1:2 or 3:4
%
% 2) Uncomment the variable you want to plot and its graph name
% plotVars_assets         = {'nrPolicy'};
% plotVarsNames_assets    = {'$i_t$'};
%
% or
%
% plotVars_assets         = {'Pii'};
% plotVarsNames_assets    = {'$\Pi_t$'};

%% Plot Figures 19 to 22 and Table 6

% Open routine setUps_polRules
%
% Uncomment the code in the section "GRAPH: Risk-free vs. Risky (Rule: target rPolicy)"
% Comment the code in other sections
% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'symbolic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate the model"
% In section "Simulation: histogram of selected variables (unconditional distribution)",
% Set graphType = 'Risk-Free vs Risky';
% Run section "Simulation: histogram of selected variables (unconditional distribution)"
%
% Uncomment the code in the section "GRAPH: Risky 4 Rules: Actual debt level"
% Comment the code in other sections
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate the model"
% In section "Simulation: histogram of selected variables (unconditional distribution)",
% Set graphType = 'Rules Comparison';
% Run section "Simulation: histogram of selected variables (unconditional distribution)"
%
% Uncomment the code in the section "GRAPH: Risky 4 Rules: 2% debt level"
% Comment the code in other sections
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate the model"
% In section "Simulation: histogram of selected variables (unconditional distribution)",
% Set graphType = 'Rules Comparison';
% Run section "Simulation: histogram of selected variables (unconditional distribution)"
%
% Uncomment the code in the section "GRAPH: Risky 4 Rules: 5% debt level"
% Comment the code in other sections
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate the model"
% In the section "Simulation: histogram of nrPolicy and pii variables (unconditional distribution)",
% Set the boolean removeNegativeNrPolicy = false;
% Run section "Simulation: histogram of nrPolicy and pii variables (unconditional distribution)"
% Set the boolean removeNegativeNrPolicy = true;
% Run section "Simulation: histogram of nrPolicy and pii variables (unconditional distribution)"

%% Plot Figures 23 to 28

% Open routine setUps_polRules
% Uncomment the code in the section "GRAPH: Risky 4 Rules: Actual debt level"
% Comment the code in other sections

% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'symbolic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate the model"
% Run section "Simulation: correlation of selected variables (unconditional distribution)"

% Repeat steps above uncommenting instead "GRAPH: Risky 4 Rules: 1% debt level"
% Repeat steps above uncommenting instead "GRAPH: Risky 4 Rules: 2% debt level"
% Repeat steps above uncommenting instead "GRAPH: Risky 4 Rules: 5% debt level"

% Run section "Simulate Correlation: default prob. and inflation"
% It shall take a long time!

% Open routine setUps_polRules
% Uncomment the code in the section "GRAPH: Risky 4 Rules: 5% debt level (fixed intercept)"
% Comment the code in other sections
% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 1;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'symbolic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% Run section "Simulate Correlation: default prob. and inflation"
% It shall take a long time!

%% Plot Figures 29 and 30

% Open routine setUps_polRules
% Uncomment the code in the section "GRAPH: Risky 4 Rules: Actual debt level"
% Comment the code in other sections

% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 2;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'automatic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% In the section "Simulate the model", set simul_order = 2;
% Run section "Simulate the model"
% Run section "Simulation: histogram of welfare variables (unconditional distribution)"

%% Plot Tables 8 to 10 and 12 to 19

% Open routine setUps_polRules
% Uncomment the code of the models rRN_inflation and polDefAdjusted_inflation in the section "GRAPH: Risky 4 Rules: Actual debt level"
% Comment the code in other sections

% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 2;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'automatic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% In the section "Simulate the model", set simul_order = 2;
% Run section "Simulate Welfare"

% Open routine setUps_polRules
% Uncomment the code of the models rRN_inflation and polDefAdjusted_inflation in the section "GRAPH: Risky 4 Rules: 1% debt level"
% Comment the code in other sections

% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 2;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'automatic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% In the section "Simulate the model", set simul_order = 2;
% Run section "Simulate Welfare"

% Open routine setUps_polRules
% Uncomment the code of the models rRN_inflation and polDefAdjusted_inflation in the section "GRAPH: Risky 4 Rules: 2% debt level"
% Comment the code in other sections

% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 2;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'automatic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% In the section "Simulate the model", set simul_order = 2;
% Run section "Simulate Welfare"

% Open routine setUps_polRules
% Uncomment the code of the models rRN_inflation and polDefAdjusted_inflation in the section "GRAPH: Risky 4 Rules: 5% debt level"
% Comment the code in other sections

% Open mainScript_model
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solveOrder  = 2;
% In the section "Parameterize, solve the model, and assign the data (Policy Rules)", set solve_derivatives_type = 'automatic';
% Run all sections until (including) section "Parameterize, solve the model, and assign the data (Policy Rules)"
% In the section "Simulate the model", set simul_order = 2;
% Run section "Simulate Welfare"