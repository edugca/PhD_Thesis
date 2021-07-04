%% Simulate correlation defPol and inflation

stateValues     = {NaN,NaN,NaN,NaN,NaN};
stateNames      = {'No solution', 'Solution is unstable', 'Solution is stable', 'Multiple solutions', 'No solution found'};

simAllMdls = {};
for iMdl = 1:length(mdlCorrVector)
    config  = riseConfigVector(iMdl);
    modelName   = riseModelNames{iMdl};
    defProcess  = riseConfigVector(iMdl).conf_defRProcessType;

    simResults = NaN(length(paramValues{1}), length(paramValues{2}));
    for iParam = 1:length(paramValues{1})
        for jParam = 1:length(paramValues{2})
            disp(['Checking stability with parameter value = ', ...
                num2str(paramValues{1}(iParam)) ' and ' num2str(paramValues{2}(jParam))]);

            paramsStructTemp = struct();
            
            if iscell(paramNames{1})
                paramsStructTemp.(paramNames{1}{1}) = paramValues{1}(iParam);
                paramsStructTemp.(paramNames{1}{2}) = paramValues{1}(iParam);
            else
                paramsStructTemp.(paramNames{1}) = paramValues{1}(iParam);
            end
            if iscell(paramNames{2})
                paramsStructTemp.(paramNames{2}{1}) = paramValues{2}(jParam);
                paramsStructTemp.(paramNames{2}{2}) = paramValues{2}(jParam);
            else
                paramsStructTemp.(paramNames{2}) = paramValues{2}(jParam);
            end

            % Override parameters
            for iParamOver = 1:length(paramNamesOverride)
                paramsStructTemp.(paramNamesOverride{iParamOver}) = paramValuesOverride{iParamOver};
            end
            mdl     = mdlCorrVector(iMdl);
            mdl     = set(mdl,'parameters', paramsStructTemp);

            % Solve model
            eigValsStr = evalc('mdl = solve(mdl);');
            disp(eigValsStr);

            % Simulate model
            simulateModel;
  
            if mdl.nsols > 0
                if is_stable_system(mdl, 'stability_algorithm', stability_algorithm)
                     % Stable solution
                     if isempty(simRecord)
                         % Simulation failed
                         simResults(iParam, jParam) = NaN;
                     else
                         simRecordCorrel = pub_CombineStructures(2,simRecord);
                         
                         % Rescale values
                         for iVar = 1:length(correlVarNames)
                            if strcmp(correlVarNames{iVar}, 'Pii')
                                correlTexVars{iVar} = '$\pi_t$ (\% annualized)';
                                y = simRecordCorrel.(correlVarNames{iVar}).data;
                                simRecordCorrel.(correlVarNames{iVar}).data = (y.^4 - 1) .* 100;
                            elseif strcmp(correlVarNames{iVar}, 'fisLim_1_2')
                                correlTexVars{iVar} = '$\mathcal{D}_t$ (\%)';
                                y = simRecordCorrel.(correlVarNames{iVar}).data;
                                simRecordCorrel.(correlVarNames{iVar}).data = y .* 100;
                            end
                         end                     
                         
                         correlSim = corr(...
                            vec(simRecordCorrel.(correlVarNames{1}).data), ...
                            vec(simRecordCorrel.(correlVarNames{2}).data) ...
                            );
                         simResults(iParam, jParam) = correlSim;
                     end
                else
                    % Non-stable solution
                    simResults(iParam, jParam) = stateValues{2};
                end
            else
                eigvals = pub_ExtractNumberFromString(eigValsStr);
                if isempty(eigvals)
                    % Regime-switching with no solution found
                    simResults(iParam, jParam) = stateValues{5};
                elseif eigvals(1) > eigvals(3)
                    % Multiple solutions
                    simResults(iParam, jParam) = stateValues{4};
                else
                    % No solution
                    simResults(iParam, jParam) = stateValues{1};
                end
            end
        end
    end
    
    simAllMdls{iMdl} = simResults;
end

%% Plot results

% Create color map
colorMapLims = [-1, 1];

% create a default color map ranging from red to light pink
colorLength = 100;
color1    = [1, 0, 0];
color2    = [1, 1, 1];
myColorMap = [linspace(color1(1),color2(1),colorLength)', linspace(color1(2),color2(2),colorLength)', linspace(color1(3),color2(3),colorLength)'];
color1    = [1, 1, 1];
color2    = [0, 0, 1];
myColorMap = [myColorMap; [linspace(color1(1),color2(1),colorLength)', linspace(color1(2),color2(2),colorLength)', linspace(color1(3),color2(3),colorLength)']];

if iscell(paramNames{1})
    str_paramNames_1 = paramNames{1}{1};
else
    str_paramNames_1 = paramNames{1};
end

if iscell(paramNames{2})
    str_paramNames_2 = paramNames{2}{1};
else
    str_paramNames_2 = paramNames{2};
end

for iMdl = 1:length(mdlCorrVector)
    
    simResults = simAllMdls{iMdl};
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
            dMatrix(iItem, 3) = simResults(ii, jj);
        end
    end

    f = figure;
    set(gca,'LooseInset',get(gca,'TightInset'));
    previousInterpreter = pub_GraphSetInterpreter('latex');
    %h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResults));
    h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResults), ...
                'MissingDataLabel', 'NA', ...
                'MissingDataColor', 'white', ...
                'FontSize', 12);
    caxis(colorMapLims);
    colormap(myColorMap);
    xlabel(paramPlotNames{2});
    ylabel(paramPlotNames{1});
    h.CellLabelFormat = '%0.2f';
    h.XDisplayLabels = num2str(str2double(h.XDisplayData),'%0.3f');
    h.YDisplayLabels = num2str(str2double(h.YDisplayData),'%0.2f');
    
    
    
    %set(gca, 'FontSize', 12);
    pub_GraphSetInterpreter(previousInterpreter);
    
%     % Add benchmark parameterization
%     ax = axes;
%     scatter(ax, paramRef{2}, paramRef{1},'LineWidth',2)
%     ax.Color = 'none'; % make invisible the background of the scatter plot

    if config.conf_hasGovernment
        simGraphName = ['withGov_correlInflationDefault_'  str_paramNames_1 '_and_'  str_paramNames_2 '_' riseConfigVector(iMdl).conf_policyRule '.png'];
    else
        simGraphName = ['withoutGov_correlInflationDefault_'  str_paramNames_1 '_and_'  str_paramNames_2 '_' riseConfigVector(iMdl).conf_policyRule '.png'];
    end
    
    set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
    exportgraphics(f, ...
    strjoin({pathImages 'Correlation' simGraphName}, filesep));

end