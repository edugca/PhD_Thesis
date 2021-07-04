%% Set-Ups

ii = 0;
clear('riseModelNames');
clear('riseSetUpNames');
clear('riseConfigVector');
riseSetUpNames = [];

%% GRAPH: ParamStability: Risky All Rules: Actual debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: fixed';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: adjusted by $\mathcal{D}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRF_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: target risky with risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: ParamStability: Risky All Rules: High debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: fixed';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: adjusted by $\mathcal{D}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRF_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: target risky with risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risk-free vs. Risky (Rule: target rPolicy)

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: adjusted by $\mathcal{D}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = true; % true, false
riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: adjusted by $\mathcal{D}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = true; % true, false    
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.01;
riseConfigVector(ii).conf_debtLevelGuess = 4.2;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: Actual debt level

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: $r^{RF}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false    
riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: adjusted by $\mathcal{D}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;
 
ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: risk-free';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRF_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: target risky with risk-free';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 1% debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: adjusted by $\mathcal{D}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRF_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: target risky with risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 2% debt level

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: $r^{RF}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false    
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.02;
riseConfigVector(ii).conf_debtLevelGuess = 4.2;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: adjusted by $\mathcal{D}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.02;
riseConfigVector(ii).conf_debtLevelGuess = 4.2;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: risk-free';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRF_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.02;
riseConfigVector(ii).conf_debtLevelGuess = 4.2;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: target risky with risk-free';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.02;
riseConfigVector(ii).conf_debtLevelGuess = 4.2;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 5% debt level

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: $r^{RF}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false    
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.05;
riseConfigVector(ii).conf_debtLevelGuess = 4.5;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: adjusted by $\mathcal{D}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.05;
riseConfigVector(ii).conf_debtLevelGuess = 4.5;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: risk-free';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRF_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.05;
riseConfigVector(ii).conf_debtLevelGuess = 4.5;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: target risky with risk-free';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.05;
riseConfigVector(ii).conf_debtLevelGuess = 4.5;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 5% debt level (fixed intercept)

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: fixed';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false    
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.05;
riseConfigVector(ii).conf_debtLevelGuess = 4.5;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% OBSERVABLES

for iConfig = 1:ii
    if riseConfigVector(iConfig).conf_estimateModel || riseConfigVector(iConfig).conf_includeObservables
    
        riseConfigVector(iConfig).conf_observablesList = [
                                ...%"obs_y", ...
                                "obs_c", ...
                                ...%"obs_pii", ...
                                ...%"obs_g", ...
                                ...%"obs_unempRate", ...
                                ...%"obs_wage", ...
                                ...%"obs_defLC", ...
                                ...%"obs_tau", ...
                                ...%"obs_debtOutput_new", ...
                                ...%"obs_nr", ...
                                ...%"obs_swap_PreDI_3m", ...
                                ...%"obs_defLC", ...
                                %"obs_y_FOCUS_20_EOQ", ...
                                %"obs_y_FOCUS_80_EOQ", ...
                                ];
        riseConfigVector(iConfig).conf_measurementErrorsList = [
                                        ...%"obs_y", ...
                                        ...%"obs_wage", ...
                                        ...%"obs_unempRate", ...
                                        ...%"obs_g", ...
                                        ...%"obs_defLC", ...
                                        ...%"obs_pii", ...
                                        %"obs_swap_PreDI_6m"
                                        %"obs_y_FOCUS_20_EOQ", ...
                                        %"obs_y_FOCUS_80_EOQ", ...
                                        %"obs_unempRate"
                                        %"obs_pii_FOCUS_Median", ...
                                        ];  
    else
    
        riseConfigVector(iConfig).conf_observablesList = [];
        riseConfigVector(iConfig).conf_measurementErrorsList = [];
        
    end
end