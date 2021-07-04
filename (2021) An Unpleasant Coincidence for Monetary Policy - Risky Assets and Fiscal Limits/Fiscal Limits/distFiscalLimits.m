%% PAPER: An Unpleasant Coincidence for Monetary Policy: Risky Assets and Fiscal Limits
%
% Author: Eduardo G. C. Amaral
% Version date: June 4th 2020
%
% This routine calculates and plots the fiscal limits

% To reestimate the fiscal limits, change the following boolean variables 
% to true in the next section of the code
% 1) reestimateFiscalLimitWithLogistic
% 2) createFiscalLimitPolicyFunctions
% 3) saveFiscalLimitPolicyFunctions

%% Simulate Fiscal Limits

%setenv('BLAS_VERBOSITY','1');
%setenv('BLAS_VERSION','c:\temp\mkl.dll');

clear all; clc;
set(0,'DefaultFigureWindowStyle','normal'); % 'normal' or 'dock'

disp('Beginning the simulation of fiscal limits');

% Final or Temporary
saveToFinalFolder = false;

% Paper
paperNumber = 2;

% Paths of the model
if paperNumber == 2
    pathMain    = edu_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) An Unpleasant Coincidence for Monetary Policy - Risky Assets and Fiscal Limits', 'C:\');
    pathData    = edu_Path('/Users/Eduardo/OneDrive/MATLAB/My Library/Database/Data/Brazil', 'C:\');
end
tablesFolder  = [pathMain filesep 'Tables' filesep 'Fiscal limits' filesep];
imagesFolder  = [pathMain filesep 'Images' filesep 'Fiscal limits' filesep];
pathSaved     = [pathMain filesep 'Saved' filesep];

% Paths of results
if saveToFinalFolder
    if paperNumber == 2
        imagesFolder = edu_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 2/Images/Fiscal limits/', 'C:\');
        tablesFolder = edu_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 2/Tables/Fiscal limits/', 'C:\');
    end
end


loadEstimParameters                 = 'All'; % 'None', 'OnlyShocks', 'All'
reestimateFiscalLimitWithLogistic   = false;
createFiscalLimitPolicyFunctions    = false;
saveFiscalLimitPolicyFunctions      = false;
loadFiscalLimitPolicyFunctions      = false;

% Define language
useLanguage = 'MatLab'; % 'C++' or 'MatLab' or 'GPU'
if strcmp(useLanguage, 'GPU')
    numVarType = 'single';
else
    numVarType = 'single';
end

%% Set up: Simulation

% Calibration name
calibrationName = 'globalsolutionModel'; % 'textbook', 'Brazil', 'globalsolutionModel'
occBinConstraint = 'FiscalLimit'; % 'None', 'FiscalLimit'
simplifyPolicyRule = false;
appendWelfareEquations = false;
stickyPrices = false;
reestimateFiscalLimitWithApproximation = false;
approximateDefaultProb = 'logistic'; % 'linear', 'logistic'
hasGovernment = true;
govBondsAreRiskFree = true;
govExpenses = 'Actual'; % 'Actual', 'High'
debtLevel = 'Actual'; % 'Actual', 'High'
govAccumDebt = true; % true, false
kVariable = 'Fixed'; % 'Fixed', 'Exogenous'
prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
policyRule = 'fixedIntercept_inflation'; % 'fixedIntercept_inflation', 'rRN_inflation', 'polDefAdjusted_inflation'
psiShockObservedWithDelay = false;
shocksTurnedOff = ["sigmaTau", "sigmaBeta"]; % 'sigmaBeta', 'sigmaA'

% Government ceiling
govCeiling = Inf; % Inf = no ceiling; simStruct.params.sigmaG = 1 std

% fixedIntercept_priceLevel,    rIntercept_priceLevel,  rRFIntercept_priceLevel
% fixedIntercept_inflation,     rIntercept_inflation,   rRFIntercept_inflation
% monPol_rule = 'fixedIntercept_inflation';

nPeriods        = 2;        % t = 1 is Steady State
nFiscalLimit    = 200;      % Number of periods at each simulation of the fiscal limit
nSimulations    = 150000;   % Number of simulations of the fiscal limit
nSimIRFs        = 1;        % Number of IRFs to be simulated
nBurnIn         = 200;      % Number of periods of shocks to be excluded
tol             = 1e-6;     % Tolerance for convergence 

% Draw fiscal limit shocks
sim_epsA    = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
sim_epsD    = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
sim_epsG    = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
sim_epsR    = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
sim_epsM    = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
sim_epsBeta = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);

% Exlude burn-in shocks
sim_epsA    = sim_epsA(:, nBurnIn+1:end);
sim_epsD    = sim_epsD(:, nBurnIn+1:end);
sim_epsG    = sim_epsG(:, nBurnIn+1:end);
sim_epsR    = sim_epsR(:, nBurnIn+1:end);
sim_epsM    = sim_epsM(:, nBurnIn+1:end);
sim_epsBeta = sim_epsBeta(:, nBurnIn+1:end);

%% List endogenous variables

% t = 1 Steady State

listFiscalLimitVars = {'sim_bbeta', 'sim_aTilde', 'sim_def', 'sim_recTilde', 'sim_g', ...
                        'sim_taxMax', 'sim_uCMax', 'sim_cMax', ...
                        'sim_yMax', 'sim_nMax', 'sim_wpMax', 'sim_tauMax', ...
                        'sim_psi', 'sim_z', 'sim_aTildeExp', 'sim_defExp', ...
                        'sim_recTildeExp', 'sim_psiExp'};


% Create vectors for the fiscal limit variables
for ii=1:length(listFiscalLimitVars)
    eval([listFiscalLimitVars{ii} ' = NaN(nSimulations,nPeriods + nFiscalLimit + 1);']);
end

% Prepare simulations
simStruct = struct();

%% Load simulations of the fiscal limit

if createFiscalLimitPolicyFunctions
    fiscalLimit_mu = NaN(11,11,11,11,11);
    fiscalLimit_std = NaN(11,11,11,11,11);
else
    load([pathSaved filesep 'estimation_fiscalLimitPolicyFunction.mat'], 'fiscalLimit_mu', 'fiscalLimit_std');
    %load('fiscalLimitPolicyFunction.mat');
end

% Load shocks
load([pathSaved filesep 'shocks.mat']);

%% Simulation at the SS (Run twice for estimation!)

disp('Simulating at the steady state');

if paperNumber == 2
    paper2_params_loadStandard;        % Parameter Values
end
params_loadValues;

%approximateFiscalLimits;

%%% Set-up 1: steady-state         
ssValues_loadStandard;  % Steady-state values
ssValues_loadValues;

simulate_fiscalLimit;

% Set paramsStruct
paramsStruct = simStruct.params;

%% Plot simulation at the SS

% Load debt to output series
sampleYearStart     = 2000;
sampleYearEnd       = 2019;

% Gross debt
% 13762	Dívida bruta do governo geral - Metodologia utilizada a partir de 2008
seriesCode = 13762;
tDebtToGDP = edu_SGS_Get(seriesCode, datetime(sampleYearStart,1,1), datetime(sampleYearEnd,1,1));
tDebtToGDP.Properties.VariableNames = {'debtToGDP'};
tDebtToGDP = retime(tDebtToGDP, 'yearly', 'lastvalue');
verticalSeries = rmmissing(tDebtToGDP);
verticalSeries = verticalSeries([datetime(2013,1,1),datetime(2016,1,1),datetime(2019,1,1)],:);

% Add IFI projections (November 17th 2020)
verticalSeries{datetime(2020,1,1), :} = 93.1;
verticalSeries{datetime(2022,1,1), :} = 97.7;

% Net debt
% 4536	Dívida líquida do governo geral (% PIB)
seriesCode = 4536;
tNetDebtToGDP = edu_SGS_Get(seriesCode, datetime(sampleYearStart,1,1), datetime(sampleYearEnd,1,1));
tNetDebtToGDP.Properties.VariableNames = {'netDebtToGDP'};
tNetDebtToGDP = retime(tNetDebtToGDP, 'yearly', 'lastvalue');
verticalSeries_netDebt = rmmissing(tNetDebtToGDP);
verticalSeries_netDebt = verticalSeries_netDebt([datetime(2019,1,1)],:);

% Plot annualized real interest rate against annual debt at the end of the period
f_r = @(x) 4 * 100 * ( ( beta.*(1-deltaBar*f_delta( 4*x, 0,0,0,0,0) ) ).^(-1) - 1 );
f_nr = @(x) 4 * 100 * ( ( beta.*(1-deltaBar*f_delta( 4*x, 0,0,0,0,0) ) ./ (1 + piiSS) ).^(-1) - 1 );
f = figure;
edu_GraphSetInterpreter('latex');
[xPlot, yPlot] = fplot(f_r, [0.0, 2.5]);
p = plot(xPlot*100, yPlot);
p.LineWidth = 3;
set(gca,'FontSize',14);
title({'Neutral real interest rate at the steady state'; ['(', edu_Num2Str(nSimulations, '.', ','), ' simulations)']}, 'interpreter', 'latex');
xlabel('steady-state debt to annual output (\%)', 'interpreter', 'latex');
ylabel('\% (annualized)', 'interpreter', 'latex');
ylim([0 30]);
xlim([0.0 150]);
yticks(0.0:5:30);
xticks(0.0:25:250);
ytickformat('%.0f');
xtickformat('%.0f');
p.Parent.XGrid = 'on';
p.Parent.YGrid = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
for iDate = 1:length(verticalSeries.Variables)
    xl = xline(verticalSeries.debtToGDP(iDate), '--r', num2str(year(verticalSeries.Properties.RowTimes(iDate))));
    xl.LineWidth = 1.5;
    xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'center';
    xl.FontWeight = 'bold';
end
for iDate = 1:length(verticalSeries_netDebt.Variables)
    xl = xline(verticalSeries_netDebt.netDebtToGDP(iDate), ':g', num2str(year(verticalSeries_netDebt.Properties.RowTimes(iDate))));
    xl.Color = [0.4660 0.6740 0.1880];
    xl.LineWidth = 2;
    xl.LabelHorizontalAlignment = 'right';
    xl.FontWeight = 'bold';
end
hold off;
%%%%%%%%%%%
simGraphName = 'fiscalLimitDist_neutralRate_ss.png';
set(gcf, 'Position',  [100, 100, 800, 400]); % resize figure
saveas(f,[imagesFolder, simGraphName]);


% Plot cdf
f = figure;
edu_GraphSetInterpreter('latex');
pdfDist = makedist('Normal', fiscalLimit_mu / (4*ySS), fiscalLimit_std / (4*ySS));
x = 0.0:0.001:2.5;
p = plot(x*100, cdf(pdfDist, x)*100);
p.LineWidth = 3;
set(gca,'FontSize',14);
title({'Fiscal limit cdf'; ['(', edu_Num2Str(nSimulations, '.', ','), ' simulations)']}, 'interpreter', 'latex');
xlabel('max debt to steady-state annual output (\%)', 'interpreter', 'latex');
ylabel('cumulative probability (\%)', 'interpreter', 'latex');
ylim([0 100]);
xlim([0.0 150]);
yticks(0.0:10:100);
xticks(0.0:25:250);
ytickformat('%.0f');
xtickformat('%.0f');
p.Parent.XGrid = 'on';
p.Parent.YGrid = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
for iDate = 1:length(verticalSeries.Variables)
    xl = xline(verticalSeries.debtToGDP(iDate), '--r', num2str(year(verticalSeries.Properties.RowTimes(iDate))));
    xl.LineWidth = 1.5;
    xl.LabelVerticalAlignment = 'middle';
    xl.LabelHorizontalAlignment = 'center';
    xl.FontWeight = 'bold';
end
for iDate = 1:length(verticalSeries_netDebt.Variables)
    xl = xline(verticalSeries_netDebt.netDebtToGDP(iDate), ':g', num2str(year(verticalSeries_netDebt.Properties.RowTimes(iDate))));
    xl.Color = [0.4660 0.6740 0.1880];
    xl.LineWidth = 2;
    xl.LabelHorizontalAlignment = 'right';
    xl.FontWeight = 'bold';
end
hold off;
%%%%%%%%%%%
simGraphName = 'fiscalLimitDist_cdf_ss.png';
set(gcf, 'Position',  [100, 100, 800, 400]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Run sensitivity analysis

% Exclude shocks
shocksExcluded  = {'epsR', 'epsD', 'epsBeta'};

% Prepate table
iTFL = 0;
tFL = array2table(NaN(11,9));
tFL.Properties.VariableNames    = {'Med_Param', 'Low_Param', 'High_Param', 'Med_Mean', 'Low_Mean', 'High_Mean', 'Med_Std', 'Low_Std', 'High_Std'};
tFL.Properties.RowNames         = { '$\overline{A}$',
                                    %'$\overline{\tilde{\omega}}$',
                                    %'$\overline{\mathcal{D}^r}$',
                                    '$\overline{G}$',
                                    '$\overline{\beta}$',
                                    '$\rho_A$',
                                    %'$\rho_R$',
                                    %'$\rho_D$',
                                    '$\rho_{GG}$',
                                    '$\rho_{GY}$',
                                    %'$\rho_{\beta}$',
                                    '$\sigma_A$',
                                    %'$\sigma_R$',
                                    %'$\sigma_D$',
                                    '$\sigma_G$',
                                    %'$\sigma_\beta$'
                                    '$\gamma_{G\Psi}$',
                                    '$\alpha_G$',
                                    '$\gamma^{NR}$'
                                  };

% Simulate fiscal limits
close all;
simulations_ss
close all;
simulations_rho
close all;
simulations_sigma
close all;
simulations_otherParams;

% Show table
%disp(tFL);

% Save table
tFL_out = tFL;
tFL_out{:,1:3} = round(tFL_out{:,1:3},3);
tFL_out{:,4:6} = round(tFL_out{:,4:6} .* 100,2);
tFL_out{:,7:9} = round(tFL_out{:,7:9} .* 100,3);

varNames = tFL_out.Properties.VariableNames;
for iVar = 1:length(varNames)
    if iVar <= 3
        tFL_out.(varNames{iVar}) = arrayfun(@(c) ThousandSep(c, '%.3f', ','),tFL_out.(varNames{iVar}),'UniformOutput',false);
    else
        tFL_out.(varNames{iVar}) = arrayfun(@(c) ThousandSep(c, '%.1f', ','),tFL_out.(varNames{iVar}),'UniformOutput',false);
    end
end

% Show formatted table
disp(tFL_out);

% In LaTex define a new column type in the beginnig of the document: \newcolumntype{Y}{>{\centering\arraybackslash}X}
colAlignment = 'l|YYY|YYY|YYY';
tabWidth = '1.0\\textwidth';
colNames = {'Medium \newline Param.', 'Low \newline Param.', 'High \newline Param.', 'Medium \newline $\mu$ \newline $(\%\, 4\overline{Y})$', 'Low \newline $\mu$ \newline $(\%\, 4\overline{Y})$', 'High \newline $\mu$ \newline $(\%\, 4\overline{Y})$', 'Medium \newline $\sigma$ \newline $(\%\, 4\overline{Y})$', 'Low \newline $\sigma$ \newline $(\%\, 4\overline{Y})$', 'High \newline $\sigma$ \newline $(\%\, 4\overline{Y})$'};
edu_Table2Latex(tFL_out, [tablesFolder, 'fiscalLimits_results.tex'], 'colAlignment', colAlignment, 'colNames', colNames, 'tabWidth', tabWidth);