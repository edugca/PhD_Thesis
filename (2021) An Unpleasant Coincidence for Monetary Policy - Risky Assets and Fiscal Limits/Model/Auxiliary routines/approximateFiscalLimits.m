%% Approximate fiscal limits by a polynomial function

% Recalculate fiscal limits for the param specification
%%%%
%%%%

% Build grid with normal processes
%nStates = 11;
if createFiscalLimitPolicyFunctions || loadFiscalLimitPolicyFunctions
    nStates = 11;
    rng_epsA        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaA,simStruct.params.rhoA,"Tauchen1986");
    rng_epsR        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaR,simStruct.params.rhoR,"Tauchen1986");
    rng_epsD        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaD,simStruct.params.rhoD,"Tauchen1986");
    rng_epsG        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaG,simStruct.params.rhoGG,"Tauchen1986"); 
    rng_epsBeta     = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaBeta,simStruct.params.rhoBeta,"Tauchen1986"); 
else
    nStates = 1;
    rng_epsA        = [0];
    rng_epsR        = [0];
    rng_epsD        = [0];
    rng_epsG        = [0];
    rng_epsBeta     = [0]; 
end
% Build grid with nodes
nodeGrid = [(1:length(rng_epsA))', ...
            (1:length(rng_epsR))', ...
            (1:length(rng_epsD))', ...
            (1:length(rng_epsG))', ...
            (1:length(rng_epsBeta))'];
gridComb  = allcomb(nodeGrid(:,1), nodeGrid(:,2), nodeGrid(:,3), nodeGrid(:,4), nodeGrid(:,5));

% Build grid with mu values
muGrid = [rng_epsA, rng_epsR, rng_epsD, rng_epsG, rng_epsBeta];
muComb = allcomb(muGrid(:,1), muGrid(:,2), muGrid(:,3), muGrid(:,4), muGrid(:,5));

% Build grid with std values
stdGrid = [rng_epsA, rng_epsR, rng_epsD, rng_epsG, rng_epsBeta];
stdComb = allcomb(stdGrid(:,1), stdGrid(:,2), stdGrid(:,3), stdGrid(:,4), stdGrid(:,5));

% Create functions
fRead_mu = @(x) fiscalLimit_mu(gridComb(x,1),gridComb(x,2),gridComb(x,3),gridComb(x,4),gridComb(x,5));
fVals_mu = arrayfun(fRead_mu, (1:size(gridComb,1))');

fRead_std = @(x) fiscalLimit_std(gridComb(x,1),gridComb(x,2),gridComb(x,3),gridComb(x,4),gridComb(x,5));
fVals_std = arrayfun(fRead_std, (1:size(gridComb,1))');

% Approximate a polynomial function for the fiscal limit mu and std
maxPower = 1;
global fitModel_mu
global fitModel_std
if nStates == 1
    %fitModel_mu     = struct();
    %fitModel_std    = struct();
    
    % Fixed point function
    fitModel_mu     = @(x1,x2, x3, x4, x5) fVals_mu;
    fitModel_std    = @(x1,x2, x3, x4, x5) fVals_std;
else
    
    % Remove NaNs of missing shocks
    idxs = find(~isnan(fVals_mu));
    muComb      = muComb(idxs, :);
    fVals_mu    = fVals_mu(idxs);
    
    idxs = find(~isnan(fVals_std));
    stdComb      = stdComb(idxs, :);
    fVals_std    = fVals_std(idxs);
    
    if isempty(muComb)
        % Fixed point function
        fitModel_mu     = @(x1,x2, x3, x4, x5) fVals_mu;
        fitModel_std    = @(x1,x2, x3, x4, x5) fVals_std;
    else 
        
        validIdx = [1,4]; % Only shocks that enter the fiscal limit!!!!
        [muComb, idxUn] = unique(muComb(:, validIdx), 'rows');
        fVals_mu    = fVals_mu(idxUn);
        [stdComb, idxUn] = unique(stdComb(:, validIdx), 'rows');
        fVals_std    = fVals_std(idxUn);
        
        sRegMu = MultiPolyRegress(muComb, fVals_mu, maxPower);
        %fitModel_mu     = @(x1,x2,x3,x4,x5) sRegMu.PolynomialExpression(x1,x4,x5);
        fMu_str = ['@(x1,x2,x3,x4,x5) ' num2str(sRegMu.Coefficients(1)) ...
                    '+' num2str(sRegMu.Coefficients(2)) '*x1 '...
                    '+' num2str(sRegMu.Coefficients(3)) '*x4 '...
                    %'+' num2str(sRegMu.Coefficients(3)) '*x5 '...
                    ];
        fitModel_mu = str2func(fMu_str);
        
        
        sRegStd = MultiPolyRegress(stdComb, fVals_std, maxPower);
        %fitModel_std     = @(x1,x2,x3,x4,x5) sRegStd.PolynomialExpression(x1,x4,x5);
        fStd_str = ['@(x1,x2,x3,x4,x5) ' num2str(sRegStd.Coefficients(1)) ...
                    '+' num2str(sRegStd.Coefficients(2)) '*x1 '...
                    '+' num2str(sRegStd.Coefficients(3)) '*x4 '...
                    %'+' num2str(sRegStd.Coefficients(4)) '*x5 '...
                    ];
        fitModel_std = str2func(fStd_str); 
        
    end
end


% Approximation of the default probability to a logistic function
if createFiscalLimitPolicyFunctions == false
    if reestimateFiscalLimitWithApproximation

        % Approximate with a logistic function
        approximateFLwithFunction();

    else

        if strcmp(approximateDefaultProb, 'logistic')
            load([pathSaved filesep 'fiscalLimit_approxLogistic.mat']);
        elseif strcmp(approximateDefaultProb, 'linear')
            load([pathSaved filesep 'fiscalLimit_approxLinear.mat']);
        end

    end
end

%fFiscalLimit_mu         = @(x1,x2,x3,x4) fitModel_mu.PolynomialExpression(x1,x2,x3,x4);
%fFiscalLimit_std        = @(x1,x2,x3,x4) fitModel_std.PolynomialExpression(x1,x2,x3,x4);
%fFiscalLimit    = @(x1,x2,x3,x4) makedist('Normal', 'mu', fFiscalLimit_mu(x1,x2,x3,x4), ...
%                                            'sigma', fFiscalLimit_std(x1,x2,x3,x4));
%fFiscalLimit_defProb    = @(bValue, x1,x2,x3,x4) normaldist(bValue, ...
%                                fFiscalLimit_mu(x1,x2,x3,x4), ...
%                                fFiscalLimit_std(x1,x2,x3,x4));

% Parameterize
if sum(muGrid(:)) + sum(stdGrid(:)) == 0 || createFiscalLimitPolicyFunctions
    simStruct.params.muFisLim_Intercept = fVals_mu;
    simStruct.params.muFisLim_aTilde    = 0;
    simStruct.params.muFisLim_g         = 0;
    simStruct.params.muFisLim_bbeta     = 0;

    simStruct.params.stdFisLim_Intercept = fVals_std;
    simStruct.params.stdFisLim_aTilde    = 0;
    simStruct.params.stdFisLim_g         = 0;
    simStruct.params.stdFisLim_bbeta     = 0;
    
    simStruct.params.probDefFisLim_Param_0      = 0;
    simStruct.params.probDefFisLim_Param_B      = 0;
    simStruct.params.probDefFisLim_Param_A      = 0;
    simStruct.params.probDefFisLim_Param_G      = 0;
    simStruct.params.probDefFisLim_Param_Beta   = 0;
else
    simStruct.params.muFisLim_Intercept = sRegMu.Coefficients(1);
    simStruct.params.muFisLim_aTilde    = sRegMu.Coefficients(2);
    simStruct.params.muFisLim_g         = sRegMu.Coefficients(3);
    %simStruct.params.muFisLim_bbeta     = sRegMu.Coefficients(4);

    simStruct.params.stdFisLim_Intercept = sRegStd.Coefficients(1);
    simStruct.params.stdFisLim_aTilde    = sRegStd.Coefficients(2);
    simStruct.params.stdFisLim_g         = sRegStd.Coefficients(3);
    %simStruct.params.stdFisLim_bbeta     = sRegStd.Coefficients(4);
    
    if strcmp(approximateDefaultProb, 'logistic')
        simStruct.params.probDefFisLim_Param_0      = mdlFL_logistic.Coefficients.Estimate(1);
        simStruct.params.probDefFisLim_Param_B      = mdlFL_logistic.Coefficients.Estimate(2);
        simStruct.params.probDefFisLim_Param_A      = mdlFL_logistic.Coefficients.Estimate(3);
        simStruct.params.probDefFisLim_Param_G      = mdlFL_logistic.Coefficients.Estimate(4);
        %simStruct.params.probDefFisLim_Param_Beta   = mdlFL_logistic.Coefficients.Estimate(7);
    elseif strcmp(approximateDefaultProb, 'linear')
        simStruct.params.probDefFisLim_Param_0      = mdlFL_linear.Coefficients.Estimate(1);
        simStruct.params.probDefFisLim_Param_B      = mdlFL_linear.Coefficients.Estimate(2);
        simStruct.params.probDefFisLim_Param_A      = mdlFL_linear.Coefficients.Estimate(3);
        simStruct.params.probDefFisLim_Param_G      = mdlFL_linear.Coefficients.Estimate(4);
        %simStruct.params.probDefFisLim_Param_Beta   = mdlFL_logistic.Coefficients.Estimate(7); 
    end
end

% % Print to file
% fMu_str = func2str(fitModel_mu);
% fMu_str = fMu_str(strfind(fMu_str, ')')+1:end);
% if strcmp(fMu_str(1), '+')
%    fMu_str = fMu_str(2:end); 
% end
% fMu_str = strrep(fMu_str, '.*', '*');
% fMu_str = strrep(fMu_str, '+-', '-');
% fMu_str = strrep(fMu_str, 'x1', 'aTilde');
% fMu_str = strrep(fMu_str, 'x2', 'recTilde');
% fMu_str = strrep(fMu_str, 'x3', 'def');
% fMu_str = strrep(fMu_str, 'x4', 'g');
% fMu_str = strrep(fMu_str, 'x5', 'bbeta');
% 
% fStd_str = func2str(fitModel_std);
% fStd_str = fStd_str(strfind(fStd_str, ')')+1:end);
% if strcmp(fStd_str(1), '+')
%    fStd_str = fStd_str(2:end); 
% end
% fStd_str = strrep(fStd_str, '.*', '*');
% fStd_str = strrep(fStd_str, '+-', '-');
% fStd_str = strrep(fStd_str, 'x1', 'aTilde');
% fStd_str = strrep(fStd_str, 'x2', 'recTilde');
% fStd_str = strrep(fStd_str, 'x3', 'def');
% fStd_str = strrep(fStd_str, 'x4', 'g');
% fStd_str = strrep(fStd_str, 'x5', 'bbeta');
% 
% fileID = fopen('rise_fiscalLimit.rs','w');
% fprintf(fileID, '%s\n', ['muFisLim = ' fMu_str ' ;']);
% fprintf(fileID, '%s\n', ['stdFisLim = ' fStd_str ' ;']);
% fclose(fileID);
% 
% fileID = fopen('rise_fiscalLimit.m','w');
% fprintf(fileID, '%s\n', ['muFisLim = ' fMu_str ' ;']);
% fprintf(fileID, '%s\n', ['stdFisLim = ' fStd_str ' ;']);
% fclose(fileID);