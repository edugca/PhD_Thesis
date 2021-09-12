%% Set-Ups

ii = 0;
clear('riseModelNames');
clear('riseSetUpNames');
clear('riseConfigVector');
riseSetUpNames = [];

%% Closed Economy; No Government; Productivity is Single Log-Normal; DefR is Null

% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'Without enterprise risk / fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 "obs_y_FOCUS_20_EOQ", ...
%                                 "obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];     

% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'Without enterprise risk / $r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 "obs_y_FOCUS_20_EOQ", ...
%                                 "obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];
                            
%% Closed Economy; No Government; Productivity is Dual Log-Normal; DefR is fixed

% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR Fixed - Fixed Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Fixed'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = '$r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR Fixed - Risky Natural Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Fixed'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false

%% Closed Economy; No Government; Productivity is Dual Log-Normal; DefR is Log-Normal

% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'With enterprise risk / fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR LogNormal - Fixed Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
                    
% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'With enterprise risk / $r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR LogNormal - Risky Natural Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
                    
%% Closed Economy; No Government; Productivity is Dual Log-Normal; DefR is MS Exogenous

% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR MS Exogenous - Fixed Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'MS_Exogenous'; % 'Null', 'Fixed', 'LogNormal', 'MS_Exogenous', 'MS_Endogenous'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 "obs_c", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 ...%"obs_y_FOCUS_20_EOQ", ...
%                                 ...%"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];
%                     
% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = '$r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR MS Exogenous - Risky Natural Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'MS_Exogenous'; % 'Null', 'Fixed', 'LogNormal', 'MS_Exogenous', 'MS_Endogenous'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 "obs_c", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 ...%"obs_y_FOCUS_20_EOQ", ...
%                                 ...%"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];

%% Closed Economy; No Government; Productivity is Dual Log-Normal; DefR is MS Endogenous

% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR MS Exogenous - Fixed Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'MS_Endogenous'; % 'Null', 'Fixed', 'LogNormal', 'MS_Exogenous', 'MS_Endogenous'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_c", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 ...%"obs_y_FOCUS_20_EOQ", ...
%                                 ...%"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];
%                     
% ii = ii + 1;
% riseModelNames{ii} = 'closedNoGovernment';
% riseSetUpNames{ii} = '$r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - No Government - Productivity Dual LogNormal - DefR MS Exogenous - Risky Natural Intercept';
% riseConfigVector(ii)                    = riseConfig;
% riseConfigVector(ii).conf_hasGovernment = false;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_policyRule    = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'MS_Endogenous'; % 'Null', 'Fixed', 'LogNormal', 'MS_Exogenous', 'MS_Endogenous'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_y_FOCUS_20_EOQ", ...
%                         "obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_c", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 ...%"obs_y_FOCUS_20_EOQ", ...
%                                 ...%"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];

%% Closed Economy; With Government; Single Regime; GovBonds are Risk-Free; Productivity is Single Log-Normal; DefR is Null

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
   
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false 
%                             
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'adjusted by $\mathcal{D}$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'target risky with RF';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
% for ii = 1:3
%     riseConfigVector(ii).conf_observablesList = [
%                             "obs_y", ...
%                             "obs_c", ...
%                             ..."obs_pii", ...
%                             "obs_g", ...
%                             "obs_unempRate", ...
%                             ...%"obs_defLC", ...
%                             ...%"obs_tau", ...
%                             ...%"obs_debtOutput_new", ...
%                             ..."obs_swap_PreDI_3m", ...
%                             %"obs_y_FOCUS_20_EOQ", ...
%                             %"obs_y_FOCUS_80_EOQ", ...
%                             ];
%     riseConfigVector(ii).conf_measurementErrorsList = [
%                                     "obs_y", ...
%                                     %"obs_g", ...
%                                     %"obs_defLC", ...
%                                     ..."obs_pii", ...
%                                     %"obs_swap_PreDI_6m"
%                                     %"obs_y_FOCUS_20_EOQ", ...
%                                     %"obs_y_FOCUS_80_EOQ", ...
%                                     "obs_unempRate",
%                                     %"obs_pii_FOCUS_Median", ...
%                                     ];  
% end

%% Closed Economy; With Government; Single Regime; Productivity is Single Log-Normal; DefR is Null

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
   
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risky: $r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false 
%                             
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'adjusted by $\mathcal{D}$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'target risky with RF';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false

% for iConfig = 1:ii
%     riseConfigVector(iConfig).conf_observablesList = [
%                             "obs_y", ...
%                             "obs_c", ...
%                             "obs_pii", ...
%                             "obs_g", ...
%                             "obs_defLC", ...
%                             ...%"obs_tau", ...
%                             ...%"obs_debtOutput_new", ...
%                             ...%"obs_swap_PreDI_3m", ...
%                             %"obs_y_FOCUS_20_EOQ", ...
%                             %"obs_y_FOCUS_80_EOQ", ...
%                             ];
%     riseConfigVector(iConfig).conf_measurementErrorsList = [
%                                     "obs_y", ...
%                                     ...%"obs_g", ...
%                                     ...%"obs_defLC", ...
%                                     ...%"obs_pii", ...
%                                     %"obs_swap_PreDI_6m"
%                                     %"obs_y_FOCUS_20_EOQ", ...
%                                     %"obs_y_FOCUS_80_EOQ", ...
%                                     %"obs_unempRate"
%                                     %"obs_pii_FOCUS_Median", ...
%                                     ];  
% end

%% GRAPH: No Debt

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Zero'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_govExpenses = 'Zero'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_govAccumDebt = false; % true, false

%% GRAPH: Estimation (OLD)

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_shocksTurnedOff = ["sigmaTau"];
% riseConfigVector(ii).conf_estimateModel = true;


% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_shocksTurnedOff = ["sigmaTau"];
% riseConfigVector(ii).conf_estimateModel = true;

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = true;

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
% riseConfigVector(ii).conf_estimateModel = true;

%% GRAPH: Estimation

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: fixed';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0;
% riseConfigVector(ii).conf_debtLevelGuess = 4;
% riseConfigVector(ii).conf_estimateModel = true;
% riseConfigVector(ii).conf_shocksTurnedOff = ["sigmaTau"];
% riseConfigVector(ii).conf_simplifyPolicyRule = false;
% riseConfigVector(ii).conf_stickyPrices = false;
% riseConfigVector(ii).conf_monopolisticCompetition = false;

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
% riseConfigVector(ii).conf_measurementErrors = true; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0;
% riseConfigVector(ii).conf_debtLevelGuess = 4;
% riseConfigVector(ii).conf_estimateModel = true;
% riseConfigVector(ii).conf_shocksTurnedOff = ["sigmaTau"];
% riseConfigVector(ii).conf_simplifyPolicyRule = false;
% riseConfigVector(ii).conf_stickyPrices = true;
% riseConfigVector(ii).conf_monopolisticCompetition = true;

%% GRAPH: Risk-free vs. Risky (Rule: target rN)

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r^{RF}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;

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
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
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
% riseConfigVector(ii).conf_measurementErrors = true; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.0125;
% riseConfigVector(ii).conf_estimateModel = false;

%% GRAPH: Risk-free vs. Risky (Rule: target rPolicy)

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risk-free: $r_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'None'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: adjusted by $\mathcal{D}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: adjusted by $\mathcal{D}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: adjusted by $\mathcal{D}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.325;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: Actual debt level

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

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: $r^{n}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rNa_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: High debt level

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

%% GRAPH: Risky 4 Rules: 1% debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: fixed';
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
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

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

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 1\%$: $r^{n}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rNa_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.01;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 2% debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: fixed';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: $r^{RF}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: adjusted by $\mathcal{D}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: risk-free';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: target risky with risk-free';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

ii = ii + 1;
riseModelNames{ii} = 'closedWithGovernment';
riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: $r^{n}_t$';
riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Natural Intercept';
riseConfigVector(ii)                        = riseConfig;
riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
riseConfigVector(ii).conf_hasGovernment     = true;
riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
riseConfigVector(ii).conf_policyRule        = 'rNa_inflation';
riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
riseConfigVector(ii).conf_measurementErrors = false; % true, false    
riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
riseConfigVector(ii).conf_defPolTarget = 0.02;
riseConfigVector(ii).conf_debtLevelGuess = 4.2;
riseConfigVector(ii).conf_estimateModel = false;
riseConfigVector(ii).conf_simplifyPolicyRule = true;

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: $r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rNaR_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 5% debt level

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
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 5\%$: $r^{n}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rNa_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.05;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Risky 4 Rules: 7% debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 7\%$: $r^{n}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rNa_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false    
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.07;
% riseConfigVector(ii).conf_debtLevelGuess = 4.8;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

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

%% GRAPH: Welfare: Raising default risk

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
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: adjusted by $\mathcal{D}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Welfare: Comparing rules at 2% default risk

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: $r^{RF}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: adjusted by $\mathcal{D}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.2;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Welfare: Comparing rules at 5% default risk

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

%% GRAPH: Welfare: Comparing rules at 7% default risk

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 7\%$: $r^{RF}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.07;
% riseConfigVector(ii).conf_debtLevelGuess = 4.7;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 7\%$: adjusted by $\mathcal{D}_t$';
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
% riseConfigVector(ii).conf_defPolTarget = 0.07;
% riseConfigVector(ii).conf_debtLevelGuess = 4.7;
% riseConfigVector(ii).conf_estimateModel = false;
% riseConfigVector(ii).conf_simplifyPolicyRule = true;

%% GRAPH: Efficient: Actual debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free';
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
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} \approx 0\%$: risk-free efficient';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free efficient';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rEF_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'Actual'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_estimateModel = false;

%% GRAPH: Efficient: High debt level

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: risk-free';
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
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: risk-free efficient';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free efficient';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rEF_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$\overline{\mathcal{D}} = 2\%$: risk-free natural';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risk-Free natural';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rNa_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = false; % true, false
% riseConfigVector(ii).conf_debtLevel = 'High'; % 'Actual', 'High', 'Zero'
% riseConfigVector(ii).conf_defPolTarget = 0.02;
% riseConfigVector(ii).conf_debtLevelGuess = 4.5;
% riseConfigVector(ii).conf_estimateModel = false;

%% Closed Economy; With Government; Productivity is Single Log-Normal; DefR is Null

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'Risky RS: $r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false    
%                             
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'adjusted by $\mathcal{D}$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'target risky with risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Single_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'Null'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% 
%                             
% for ii = 1:4
%     riseConfigVector(ii).conf_observablesList = [
%                             "obs_y", ...
%                             %"obs_c", ...
%                             %"obs_pii", ...
%                             %"obs_g", ...
%                             %"obs_defLC", ...
%                             %"obs_tau", ...
%                             %"obs_debtOutput_new", ...
%                             %"obs_swap_PreDI_3m", ...
%                             %"obs_y_FOCUS_20_EOQ", ...
%                             %"obs_y_FOCUS_80_EOQ", ...
%                             ];
%     riseConfigVector(ii).conf_measurementErrorsList = [
%                                     "obs_y", ...
%                                     ...%"obs_g", ...
%                                     %"obs_defLC", ...
%                                     %"obs_pii", ...
%                                     %"obs_swap_PreDI_6m"
%                                     %"obs_y_FOCUS_20_EOQ", ...
%                                     %"obs_y_FOCUS_80_EOQ", ...
%                                     %"obs_unempRate"
%                                     %"obs_pii_FOCUS_Median", ...
%                                     ];  
% end

%% Closed Economy; With Government; Productivity is Dual Log-Normal; DefR is Log-Normal; Gov Bonds are Risk-Free

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         %"obs_y", ...
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];     
% 
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         %"obs_y", ...
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];     
%                             
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'adjusted by $\mathcal{D}$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = true; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         %"obs_y", ...
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];
                            
%% Closed Economy; With Government; Productivity is Dual Log-Normal; DefR is Log-Normal; Gov Bonds are Risky

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'fixed intercept';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Fixed Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'fixedIntercept_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         "obs_y", ...
%                         "obs_c", ...
%                         "obs_swap_PreDI_3m", ...
%                         "obs_pii", ...
%                         "obs_defLC", ...
%                         %"obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 "obs_y", ... 
%                                 "obs_swap_PreDI_3m", ...
%                                 "obs_pii", ...
%                                 "obs_defLC", ...
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];     

% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = '$r^{rn}_t$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - Risky Natural Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRN_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         %"obs_y", ...
%                         "obs_c", ...
%                         %"obs_swap_PreDI_3m", ...
%                         %"obs_pii", ...
%                         %"obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];     
%                             
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'adjusted by $\mathcal{D}$';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'polDefAdjusted_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         %"obs_y", ...
%                         "obs_c", ...
%                         %"obs_swap_PreDI_3m", ...
%                         %"obs_pii", ...
%                         %"obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];
%                             
% ii = ii + 1;
% riseModelNames{ii} = 'closedWithGovernment';
% riseSetUpNames{ii} = 'target risky with risk-free';
% riseFullSetUpNames{ii} = 'Closed Economy - With Government - Productivity Single LogNormal - DefR Null - PolDef Adjusted Intercept';
% riseConfigVector(ii)                        = riseConfig;
% riseConfigVector(ii).conf_occBinConstraint  = 'FiscalLimit'; % 'None', 'FiscalLimit'
% riseConfigVector(ii).conf_hasGovernment     = true;
% riseConfigVector(ii).conf_govBondsAreRiskFree = false; % true, false
% riseConfigVector(ii).conf_policyRule        = 'rRNwithRiskFree_inflation';
% riseConfigVector(ii).conf_prodProcessType = 'Dual_LogNormal'; % 'Single_LogNormal', 'Dual_LogNormal', 'MarkovSwitching'
% riseConfigVector(ii).conf_defRProcessType = 'LogNormal'; % 'Null', 'Fixed', 'LogNormal'
% riseConfigVector(ii).conf_measurementErrors = true; % true, false
% riseConfigVector(ii).conf_observablesList = [
%                         %"obs_y", ...
%                         "obs_c", ...
%                         %"obs_swap_PreDI_3m", ...
%                         %"obs_pii", ...
%                         %"obs_g", ...
%                         %"obs_y_FOCUS_20_EOQ", ...
%                         %"obs_y_FOCUS_80_EOQ", ...
%                         ];
% riseConfigVector(ii).conf_measurementErrorsList = [  
%                                 ...%"obs_y", ...                              
%                                 ...%"obs_pii", ...
%                                 ...%"obs_swap_PreDI_6m"
%                                 %"obs_y_FOCUS_20_EOQ", ...
%                                 %"obs_y_FOCUS_80_EOQ", ...
%                                 %"obs_unempRate"
%                                 %"obs_pii_FOCUS_Median", ...
%                                 ];   

%% OBSERVABLES

for iConfig = 1:ii
    if riseConfigVector(iConfig).conf_estimateModel || riseConfigVector(iConfig).conf_includeObservables
    
        riseConfigVector(iConfig).conf_observablesList = [
                                "obs_y", ...
                                "obs_c", ...
                                "obs_pii", ...
                                "obs_g", ...
                                ..."obs_unempRate", ...
                                ..."obs_wage", ...
                                ..."obs_defLC", ...
                                ..."obs_tau", ...
                                ..."obs_debtOutput_new", ...
                                ..."obs_nr", ...
                                "obs_swap_PreDI_3m", ...
                                ..."obs_defLC", ...
                                ..."obs_y_FOCUS_20_EOQ", ...
                                ..."obs_y_FOCUS_80_EOQ", ...
                                ];
        riseConfigVector(iConfig).conf_measurementErrorsList = [
                                        "obs_y", ...
                                        ...%"obs_c", ...
                                        ...%"obs_wage", ...
                                        ...%"obs_unempRate", ...
                                        ...%"obs_g", ...
                                        ..."obs_defLC", ...
                                        ...%"obs_pii", ...
                                        ..."obs_swap_PreDI_6m"
                                        ..."obs_y_FOCUS_20_EOQ", ...
                                        ..."obs_y_FOCUS_80_EOQ", ...
                                        ..."obs_unempRate"
                                        ..."obs_pii_FOCUS_Median", ...
                                        ];  
    else
    
        riseConfigVector(iConfig).conf_observablesList = [];
        riseConfigVector(iConfig).conf_measurementErrorsList = [];
        
    end
end