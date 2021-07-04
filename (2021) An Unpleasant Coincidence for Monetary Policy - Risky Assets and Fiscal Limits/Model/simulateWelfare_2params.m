
simResults   = NaN(length(paramValues{1}), length(paramValues{2}));
simResultsR  = NaN(length(paramValues{1}), length(paramValues{2}));
simResultsNR = NaN(length(paramValues{1}), length(paramValues{2}));

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
        mdl     = mdlSimVector(iMdl);
        mdl     = set(mdl,'parameters', paramsStructTemp);

        % Solve model
        eigValsStr = evalc('mdl = solve(mdl);');
        disp(eigValsStr);

        if mdl.nsols > 0
            if is_stable_system(mdl, 'stability_algorithm', stability_algorithm)
                 % Simulate model
                 simulateModel;

                 % Check if there were more than 1 simulation
                 if size(simRecord, 1) > 1
                    data_welfareR = [simRecord(:).welfareR];
                    data_welfareNR = [simRecord(:).welfareNR];
                 else
                     data_welfareR = [simRecord.welfareR];
                     data_welfareNR = [simRecord.welfareNR];
                 end
                 vec_welfareR = vec(data_welfareR.data);
                 vec_welfareNR = vec(data_welfareNR.data);

                 % Stable solution
                 % Sum welfare of Ricardian and non-Ricardian agents
                 if isempty(simRecord)
                     % Simulation failed
                     simResults(iParam, jParam)   = NaN;
                     simResultsR(iParam, jParam)  = NaN;
                     simResultsNR(iParam, jParam) = NaN;
                 else
                    simResults(iParam, jParam) = ...
                        mean(...
                        (1-paramsStruct.fracNR)*vec_welfareR ...
                        + paramsStruct.fracNR*vec_welfareNR ...
                        );
                    simResultsR(iParam, jParam) = ...
                        mean(vec_welfareR);
                    simResultsNR(iParam, jParam) = ...
                        mean(vec_welfareNR);
                 end
            else
                % Non-stable solution
                simResults(iParam, jParam)   = stateValues{2};
                simResultsR(iParam, jParam)  = stateValues{2};
                simResultsNR(iParam, jParam) = stateValues{2};
            end
        else
            eigvals = edu_ExtractNumberFromString(eigValsStr);
            if isempty(eigvals)
                % Regime-switching with no solution found
                simResults(iParam, jParam)   = stateValues{5};
                simResultsR(iParam, jParam)  = stateValues{5};
                simResultsNR(iParam, jParam) = stateValues{5};
            elseif ~isempty(eigvals) && eigvals(1) > eigvals(3)
                % Multiple solutions
                simResults(iParam, jParam)   = stateValues{4};
                simResultsR(iParam, jParam)  = stateValues{4};
                simResultsNR(iParam, jParam) = stateValues{4};
            else
                % No solution
                simResults(iParam, jParam)   = stateValues{1};
                simResultsR(iParam, jParam)  = stateValues{1};
                simResultsNR(iParam, jParam) = stateValues{1};
            end
        end
    end
end

simAllMdls{iMdl}    = simResults;
simAllMdls_R{iMdl}  = simResultsR;
simAllMdls_NR{iMdl} = simResultsNR;