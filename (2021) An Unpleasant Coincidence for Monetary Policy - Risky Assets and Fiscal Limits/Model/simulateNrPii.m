%% Simulate Interest Rates and Inflation

close all;

% Initiate randomizer
rng(1900)

stateValues     = {NaN,NaN,NaN,NaN,NaN};
stateNames      = {'No solution', 'Solution is unstable', 'Solution is stable', 'Multiple solutions', 'RS with no solution found'};

simAllMdls = {};
%simAllMdls_R = {};
%simAllMdls_NR = {};
for iMdl = 1:length(mdlSimVector)
    config  = riseConfigVector(iMdl);
    modelName   = riseModelNames{iMdl};
    defProcess  = riseConfigVector(iMdl).conf_defRProcessType;

    if length(paramValues) == 2
        simulateNrPii_2params;
    elseif length(paramValues) == 3
        simulateNrPii_3params;
    end
end

%% Table

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

% Model used as reference against which others will be compared
refMdl = 2;

% Select plot params
iParams = find(ismember(paramValues{1}, 1:0.0250:3.0));
jParams = find(ismember(paramValues{2}, 0.0:0.15:0.9));
kParams = find(ismember(paramValues{3}, 0.0));

rowNames = string(strcat([paramTexPlotNames{1} '='], num2str(paramValues{1}','%.2f')));
varNames = string(strcat([paramTexPlotNames{2} '='], num2str(paramValues{2}','%.3f')));

tabWidth = '0.9\\textwidth';
colAlignment = repmat('r', 1, length(varNames));
nRound = 1;

for iMdl = 1:length(mdlSimVector)
    
    for kParam = 1:length(kParams)
    
        t_nrPolicy = array2table(round(simAllMdls{iMdl}.nrPolicy,nRound), 'VariableNames', varNames, 'RowNames', rowNames);
        t_rPolicy = array2table(round(simAllMdls{iMdl}.rPolicy,nRound), 'VariableNames', varNames, 'RowNames', rowNames);
        t_Pii = array2table(round(simAllMdls{iMdl}.Pii,nRound), 'VariableNames', varNames, 'RowNames', rowNames);
        t_rNa = array2table(round(simAllMdls{iMdl}.rNa,nRound), 'VariableNames', varNames, 'RowNames', rowNames);
        t_probDefFisLim = array2table(round(simAllMdls{iMdl}.probDefFisLim,nRound), 'VariableNames', varNames, 'RowNames', rowNames);

        simGraphName = ['Rates and Inflation - XXXXXX - ' mdlSimVector(iMdl).user_data.conf_policyRule ' - ' ...
            mdlSimVector(iMdl).user_data.conf_debtLevel ' - ' num2str(mdlSimVector(iMdl).user_data.conf_defPolTarget) ...
            ' - ' paramNames{3} '=' num2str(paramValues{3}(kParam)) ...
            pub_IIF(isempty(paramNamesOverride), '', @() [' - ' paramNamesOverride{1} '=' num2str(paramValuesOverride{1})]) ...
            '.tex'];


        % nrPolicy
        simGraphName_nrPolicy = replace(simGraphName, 'XXXXXX', 'nrPolicy');
        filePath = strjoin({pathTables 'Simulation' 'Rates and Inflation' title_stickyPrices simGraphName_nrPolicy}, filesep);
        pub_Table2Latex(t_nrPolicy(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);

         % rPolicy
        simGraphName_rPolicy = replace(simGraphName, 'XXXXXX', 'rPolicy');
        filePath = strjoin({pathTables 'Simulation' 'Rates and Inflation' title_stickyPrices simGraphName_rPolicy}, filesep);
        pub_Table2Latex(t_rPolicy(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);

         % Pii
        simGraphName_Pii = replace(simGraphName, 'XXXXXX', 'Pii');
        filePath = strjoin({pathTables 'Simulation' 'Rates and Inflation' title_stickyPrices simGraphName_Pii}, filesep);
        pub_Table2Latex(t_Pii(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);
    
        % rNa
        simGraphName_rNa = replace(simGraphName, 'XXXXXX', 'rNa');
        filePath = strjoin({pathTables 'Simulation' 'Rates and Inflation' title_stickyPrices simGraphName_rNa}, filesep);
        pub_Table2Latex(t_rNa(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);
        
        % rNa
        simGraphName_probDefFisLim = replace(simGraphName, 'XXXXXX', 'probDefFisLim');
        filePath = strjoin({pathTables 'Simulation' 'Rates and Inflation' title_stickyPrices simGraphName_probDefFisLim}, filesep);
        pub_Table2Latex(t_probDefFisLim(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);
        
    end
    
end


%% Welfare extra (NOT BEING USED)
%simulateWelfareExtra;
%%%%%
