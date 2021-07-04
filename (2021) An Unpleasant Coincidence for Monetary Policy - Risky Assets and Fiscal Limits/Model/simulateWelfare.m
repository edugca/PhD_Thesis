%% Simulate welfare

close all;

% Initiate randomizer
rng(1900)

stateValues     = {NaN,NaN,NaN,NaN,NaN};
stateNames      = {'No solution', 'Solution is unstable', 'Solution is stable', 'Multiple solutions', 'RS with no solution found'};

simAllMdls = {};
simAllMdls_R = {};
simAllMdls_NR = {};
for iMdl = 1:length(mdlSimVector)
    config  = riseConfigVector(iMdl);
    modelName   = riseModelNames{iMdl};
    defProcess  = riseConfigVector(iMdl).conf_defRProcessType;

    if length(paramValues) == 2
        simulateWelfare_2params;
    elseif length(paramValues) == 3
        simulateWelfare_3params;
    end
end

%% Table welfare

if stickyPrices
   title_stickyPrices = 'Sticky Prices'; 
else
   title_stickyPrices = 'Flexible Prices'; 
end

% Model used as reference against which others will be compared
refMdl = 2;

% Select plot params
iParams = find(ismember(paramValues{1}, 1:0.0250:3.0));
jParams = find(ismember(paramValues{2}, 0.10:0.0250:0.30));

rowNames = string(strcat([paramTexPlotNames{1} '='], num2str(paramValues{1}','%.2f')));
varNames = string(strcat([paramTexPlotNames{2} '='], num2str(paramValues{2}','%.3f')));

tabWidth = '0.75\\textwidth';
colAlignment = repmat('r', 1, length(varNames));
nRound = 0;

for iMdl = 1:length(mdlSimVector)

    t_Aggregate = array2table(round(simAllMdls{iMdl},nRound), 'VariableNames', varNames, 'RowNames', rowNames);
    t_Ricardian = array2table(round(simAllMdls_R{iMdl},nRound), 'VariableNames', varNames, 'RowNames', rowNames);
    t_NonRicardian = array2table(round(simAllMdls_NR{iMdl},nRound), 'VariableNames', varNames, 'RowNames', rowNames);

    % Aggregate
    simGraphName = ['Welfare - Aggregate -' mdlSimVector(iMdl).user_data.conf_policyRule ' - ' ...
        mdlSimVector(iMdl).user_data.conf_debtLevel ' - ' num2str(mdlSimVector(iMdl).user_data.conf_defPolTarget) '.tex'];
    filePath = strjoin({pathTables 'Welfare' title_stickyPrices simGraphName}, filesep);
    pub_Table2Latex(t_Aggregate(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);

    % Ricardian
    simGraphName = ['Welfare - Ricardian -' mdlSimVector(iMdl).user_data.conf_policyRule ' - ' ...
        mdlSimVector(iMdl).user_data.conf_debtLevel ' - ' num2str(mdlSimVector(iMdl).user_data.conf_defPolTarget) '.tex'];
    filePath = strjoin({pathTables 'Welfare' title_stickyPrices simGraphName}, filesep);
    pub_Table2Latex(t_Ricardian(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);

    % NonRicardian
    simGraphName = ['Welfare - NonRicardian -' mdlSimVector(iMdl).user_data.conf_policyRule ' - ' ...
        mdlSimVector(iMdl).user_data.conf_debtLevel ' - ' num2str(mdlSimVector(iMdl).user_data.conf_defPolTarget) '.tex'];
    filePath = strjoin({pathTables 'Welfare' title_stickyPrices simGraphName}, filesep);
    pub_Table2Latex(t_NonRicardian(iParams,jParams), filePath, 'tabWidth', tabWidth, 'colAlignment', colAlignment);

end


%% Welfare extra (NOT BEING USED)
%simulateWelfareExtra;
%%%%%
