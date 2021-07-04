%% Simulate

% t = 1 Steady State

% Check whether it is to create policyFunction vectors
if createFiscalLimitPolicyFunctions
    startPoint = 'SS';
    
    % Draw shocks
    epsA = zeros(nSimIRFs,nPeriods + nFiscalLimit + 1, numVarType);
    epsD = zeros(nSimIRFs,nPeriods + nFiscalLimit + 1, numVarType);
    epsG = zeros(nSimIRFs,nPeriods + nFiscalLimit + 1, numVarType);
    epsR = zeros(nSimIRFs,nPeriods + nFiscalLimit + 1, numVarType);
    epsM = zeros(nSimIRFs,nPeriods + nFiscalLimit + 1, numVarType);
    epsBeta = zeros(nSimIRFs,nPeriods + nFiscalLimit + 1, numVarType);

    % Draw fiscal limit shocks
    sim_epsA = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
    sim_epsD = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
    sim_epsG = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
    sim_epsR = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
    sim_epsM = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);
    sim_epsBeta = norminv(rand(nSimulations,nPeriods + nFiscalLimit + nBurnIn + 1, numVarType),0,1);

    % Exlude burn-in shocks
    sim_epsA = sim_epsA(:, nBurnIn+1:end);
    sim_epsD = sim_epsD(:, nBurnIn+1:end);
    sim_epsG = sim_epsG(:, nBurnIn+1:end);
    sim_epsR = sim_epsR(:, nBurnIn+1:end);
    sim_epsM = sim_epsM(:, nBurnIn+1:end);
    sim_epsBeta = sim_epsBeta(:, nBurnIn+1:end);
    
    %epsA(:,:) = 0;
    %epsR(:,:) = 0;
    %epsD(:,:) = 0;
    %epsG(:,:) = 0;

    if paperNumber == 2
        paper2_params_loadStandard;        % Parameter Values
    elseif paperNumber == 3
        paper3_params_loadStandard;        % Parameter Values
    end
    params_loadValues;

    %%% Set-up 1: steady-state         
    ssValues_loadStandard;  % Steady-state values
    ssValues_loadValues;

    nStates = 11;
    if strcmp(prodProcessType, 'Single_LogNormal')
        rng_epsA    = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaA, simStruct.params.rhoA,"Tauchen1986");
        rng_epsR    = [0];
        rng_epsD    = [0];
        rng_epsG    = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaG, simStruct.params.rhoGG,"Tauchen1986");
        rng_epsM    = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaM, simStruct.params.rhoM,"Tauchen1986"); 
        rng_epsBeta = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaBeta, simStruct.params.rhoBeta,"Tauchen1986"); 
    else
        rng_epsA        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaA, simStruct.params.rhoA,"Tauchen1986");
        rng_epsR        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaR, simStruct.params.rhoR,"Tauchen1986");
        rng_epsD        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaD, simStruct.params.rhoD,"Tauchen1986");
        rng_epsG        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaG, simStruct.params.rhoGG,"Tauchen1986");
        rng_epsM        = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaM, simStruct.params.rhoM,"Tauchen1986");
        rng_epsBeta     = edu_arToMarkov(nStates,3,'Normal',0,simStruct.params.sigmaBeta, simStruct.params.rhoBeta,"Tauchen1986");
    end
    
    % Turn off shocks
    for iShock = 1:length(shocksTurnedOff)
        eval(['rng_eps' shocksTurnedOff{iShock}(6:end) ' = [0]' ';']);
    end
    
else
    nStates = 1;
    rng_epsA    = zeros(1, 1, numVarType);
    rng_epsR    = zeros(1, 1, numVarType);
    rng_epsD    = zeros(1, 1, numVarType);
    rng_epsG    = zeros(1, 1, numVarType);
    rng_epsM    = zeros(1, 1, numVarType);
    rng_epsBeta = zeros(1, 1, numVarType);
end

% List of simulation variables
listFiscalLimitVars = {'sim_aTilde', 'sim_def', 'sim_recTilde', 'sim_g', ...
                        'sim_taxMax', 'sim_uCMax', 'sim_cMax', ...
                        'sim_yMax', 'sim_nMax', 'sim_wpMax', 'sim_tauMax', ...
                        'sim_Psi', 'sim_z', 'sim_tLS', 'sim_aTildeExp', 'sim_defExp', ...
                        'sim_recTildeExp', 'sim_PsiExp', 'sim_pii'};

%% Pre-allocate for parallelization (keeps simulation path)
% sim_aTilde      = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_def         = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_recTilde    = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_g           = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_taxMax      = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_uCMax       = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_cMax        = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_yMax        = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_nMax        = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% %sim_wpMax       = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_tauMax      = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_psi         = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_z           = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_tLS         = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% 
% sim_aTildeExp   = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_defExp      = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_recTildeExp = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_psiExp      = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% sim_gExp        = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));
% 
% %sim_pii         = complex(NaN(nSimulations,nPeriods + nFiscalLimit + 1, numVarType));


%% Pre-allocate for parallelization (does not keep simulation path)

if ~strcmp(useLanguage, 'MatLab')
    rng_epsA    = complex(rng_epsA);
    rng_epsR    = complex(rng_epsR);
    rng_epsD    = complex(rng_epsD);
    rng_epsG    = complex(rng_epsG);
    rng_epsM    = complex(rng_epsM);
    rng_epsBeta = complex(rng_epsBeta);
    
    sim_epsA    = complex(sim_epsA);
    sim_epsR    = complex(sim_epsR);
    sim_epsD    = complex(sim_epsD);
    sim_epsG    = complex(sim_epsG);
    sim_epsM    = complex(sim_epsM);
    sim_epsBeta = complex(sim_epsBeta);
    
    sim_aTilde      = complex(NaN(nSimulations,2, numVarType));
    sim_def         = complex(NaN(nSimulations,2, numVarType));
    sim_recTilde    = complex(NaN(nSimulations,2, numVarType));
    sim_g           = complex(NaN(nSimulations,2, numVarType));
    sim_taxMax      = complex(NaN(nSimulations,2, numVarType));
    sim_uCMax       = complex(NaN(nSimulations,2, numVarType));
    sim_cMax        = complex(NaN(nSimulations,2, numVarType));
    sim_yMax        = complex(NaN(nSimulations,2, numVarType));
    sim_nMax        = complex(NaN(nSimulations,2, numVarType));
    %sim_wpMax       = complex(NaN(nSimulations,2, numVarType));
    sim_tauMax      = complex(NaN(nSimulations,2, numVarType));
    sim_psi         = complex(NaN(nSimulations,2, numVarType));
    sim_z           = complex(NaN(nSimulations,2, numVarType));
    sim_tLS         = complex(NaN(nSimulations,2, numVarType));

    sim_aTildeExp   = complex(NaN(nSimulations,2, numVarType));
    sim_defExp      = complex(NaN(nSimulations,2, numVarType));
    sim_recTildeExp = complex(NaN(nSimulations,2, numVarType));
    sim_psiExp      = complex(NaN(nSimulations,2, numVarType));
    sim_gExp        = complex(NaN(nSimulations,2, numVarType));

    %sim_pii         = complex(NaN(nSimulations,2, numVarType));

end



%%

% Clean vectors for the fiscal limit variables
% and attribute fist value to steady-state value
for ii=1:length(listFiscalLimitVars)
    eval([listFiscalLimitVars{ii} '(:,:) = NaN;']);
    eval([listFiscalLimitVars{ii} '(:,1) = ' listFiscalLimitVars{ii}(5:end) 'SS;']);
end

% Steady-state
sim_epsA(:,1)       = 0;
sim_epsR(:,1)       = 0;
sim_epsD(:,1)       = 0;
sim_epsG(:,1)       = 0;
sim_epsM(:,1)       = 0;
sim_epsBeta(:,1)    = 0;

% CONSTANT SERIES
sim_z(:,:)          = zSS;
sim_tLS(:,:)        = tLSSS;
sim_tauMax(:,:)     = tauMaxSS; %chi/(1+chi)

% Call Fiscal Limit Simulation
caller_simulateFiscalLimit_CCode;

% Save fiscal limit policy function
if saveFiscalLimitPolicyFunctions
    
    %%%%% Extrapolate over NaNs
    if sum(isnan(fiscalLimit_mu(:)))>0
        for iR = 1:size(fiscalLimit_mu,2)
            fiscalLimit_mu(:,iR,:,:,:) = fiscalLimit_mu(:,1,:,:,:);
            fiscalLimit_std(:,iR,:,:,:) = fiscalLimit_std(:,1,:,:,:);
        end
        for iD = 1:size(fiscalLimit_mu,2)
            fiscalLimit_mu(:,:,iD,:,:) = fiscalLimit_mu(:,:,1,:,:);
            fiscalLimit_std(:,:,iD,:,:) = fiscalLimit_std(:,:,1,:,:);
        end
    end
    
    if strcmp(loadEstimParameters, 'All')
        save([pathSaved filesep 'estimation_fiscalLimitPolicyFunction.mat'], 'fiscalLimit_mu', 'fiscalLimit_std');
    elseif strcmp(loadEstimParameters, 'OnlyShocks')
        save([pathSaved filesep 'estimOnlyShocks_fiscalLimitPolicyFunction.mat'], 'fiscalLimit_mu', 'fiscalLimit_std');
    elseif strcmp(loadEstimParameters, 'None')
        save([pathSaved filesep 'calibration_fiscalLimitPolicyFunction.mat'], 'fiscalLimit_mu', 'fiscalLimit_std');
    end

end