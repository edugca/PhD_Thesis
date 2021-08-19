
simResults   = struct();
simResults.nrPolicy     = NaN(length(paramValues{1}), length(paramValues{2}));
simResults.rPolicy      = NaN(length(paramValues{1}), length(paramValues{2}));
simResults.Pii          = NaN(length(paramValues{1}), length(paramValues{2}));

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
                    data_nrPolicy   = [simRecord(:).nrPolicy];
                    data_rPolicy    = [simRecord(:).rPolicy];
                    data_Pii        = [simRecord(:).Pii];
                 else
                     data_nrPolicy  = [simRecord.nrPolicy];
                     data_rPolicy   = [simRecord.rPolicy];
                     data_Pii       = [simRecord.Pii];
                 end
                 vec_nrPolicy   = vec(data_nrPolicy.data);
                 vec_rPolicy    = vec(data_rPolicy.data);
                 vec_Pii        = vec(data_Pii.data);

                 % Stable solution
                 % Sum welfare of Ricardian and non-Ricardian agents
                 if isempty(simRecord)
                     % Simulation failed
                     simResults.nrPolicy(iParam, jParam)    = NaN;
                     simResults.rPolicy(iParam, jParam)     = NaN;
                     simResults.Pii(iParam, jParam)         = NaN;
                 else
                    simResults.nrPolicy(iParam, jParam) = mean( ( (1 + vec_nrPolicy) .^ 4 - 1 ) .* 100 );
                    simResults.rPolicy(iParam, jParam) = mean( ( (1 + vec_rPolicy) .^ 4 - 1 ) .* 100 );
                    simResults.Pii(iParam, jParam) = mean( ( vec_Pii .^ 4 - 1 ) .* 100 );
                 end
            else
                % Non-stable solution
                simResults.nrPolicy(iParam, jParam)     = stateValues{2};
                simResults.rPolicy(iParam, jParam)      = stateValues{2};
                simResults.Pii(iParam, jParam)          = stateValues{2};
            end
        else
            eigvals = edu_ExtractNumberFromString(eigValsStr);
            if isempty(eigvals)
                % Regime-switching with no solution found
                simResults(iParam, jParam).nrPolicy   = stateValues{5};
                simResults(iParam, jParam).rPolicy   = stateValues{5};
                simResults(iParam, jParam).Pii   = stateValues{5};
            elseif ~isempty(eigvals) && eigvals(1) > eigvals(3)
                % Multiple solutions
                simResults(iParam, jParam).nrPolicy   = stateValues{4};
                simResults(iParam, jParam).rPolicy    = stateValues{4};
                simResults(iParam, jParam).Pii        = stateValues{4};
            else
                % No solution
                simResults(iParam, jParam).nrPolicy     = stateValues{1};
                simResults(iParam, jParam).rPolicy      = stateValues{1};
                simResults(iParam, jParam).Pii          = stateValues{1};
            end
        end
    end
end

simAllMdls{iMdl}    = simResults;
