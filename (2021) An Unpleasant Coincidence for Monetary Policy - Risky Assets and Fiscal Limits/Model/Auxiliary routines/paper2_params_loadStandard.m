%% Parameter Values

%% Basic parameters
simStruct.params.beta       = 0.989;    % Castro (2015)
simStruct.params.sigma      = 1.3;      % Castro (2015)
simStruct.params.chi        = 1;        % Castro (2015)
simStruct.params.nBar       = 1/3;      % Traditional in the literature
simStruct.params.yBar       = 1;        % Just a normalization
simStruct.params.ySS        = simStruct.params.yBar;
simStruct.params.kBar       = 18.0 * simStruct.params.yBar; % PennWorld table
simStruct.params.ratioRecA  = 0.8;      % Just a normalization
simStruct.params.defBar     = 0.5;      % Just a normalization
simStruct.params.pLevelBar  = 1;        % Just a normalization
simStruct.params.fracNR     = 0.40;     % Castro (2015)
simStruct.params.alphaG     = 0;        % -1.51 Fève and Sahuc (2017); (POOR) Estimation: 0.437
simStruct.params.gammaGPSI  = 0.2;      % Estimation: 0.510

% Override with Galí Textbook, Chapter 3
if strcmp(calibrationName, 'textbook')
    
    simStruct.params.beta      = 0.99;
    simStruct.params.sigma     = 1.0;
    simStruct.params.chi       = 1.0;
    simStruct.params.nBar      = 1/3;
    simStruct.params.yBar      = 1;
    simStruct.params.kBar      = 1;
    
    simStruct.params.ratioRecA  = 0.8;
    simStruct.params.defBar     = 0.05;
    
elseif  strcmp(calibrationName, 'Brazil')
    
    simStruct.params.beta      = 0.989; % Castro (2015)
    simStruct.params.sigma     = 1.3;   % Castro (2015)
    simStruct.params.chi       = 1.0;   % Castro (2015)
    simStruct.params.nBar      = 1/3;   % Traditional in the literature
    simStruct.params.yBar      = 1;     % Just a normalization
    simStruct.params.kBar      = 1;     % Just a normalization
    
    simStruct.params.ratioRecA  = 0.8;  % Just stylized
    simStruct.params.defBar     = 0.05; % Just stylized
end

if strcmp(defRProcessType, 'Null')
     simStruct.params.defBar     = 0;
     simStruct.params.ratioRecA  = 0;
end

%% STICKY PRICES

simStruct.params.thetaElast = 11;
simStruct.params.labSub     = 0; % labor subsidy: i.e. 1/thetaElast

if stickyPrices
    % Sticky prices
    simStruct.params.phiC       = 100;
else
    % Flexible prices
    simStruct.params.phiC       = 0;
end

%% SHOCK PARAMETERS
simStruct.params.rhoA       = 0.964; % Estimation
simStruct.params.rhoR       = 0.964; % Arbitary
simStruct.params.rhoGG      = 0.951; % Estimation
simStruct.params.rhoGY      = 0.036; % Estimation
simStruct.params.rhoTau     = 0.862; % Estimation (non-Bayesian)
simStruct.params.rhoD       = 0.5;  % Arbitrary
simStruct.params.rhoM       = 0.5;  % Arbitrary
simStruct.params.rhoBeta    = 0.666; % Estimation
simStruct.params.rhoPolDef  = 0.5;  % Arbitrary

simStruct.params.rhoAG       = 0.0; % Arbitrary

simStruct.params.sigmaA      = 0.007;  % Estimation
simStruct.params.sigmaR      = 0.007;  % Arbitrary
simStruct.params.sigmaG      = 0.012;  % Estimation
simStruct.params.sigmaTau    = 0.00938; % Estimation (non-Bayesian)
simStruct.params.sigmaD      = 0.0100; % Arbitrary
simStruct.params.sigmaM      = 0.003;  % Estimation
simStruct.params.sigmaPsi    = 0.0100; % Arbitrary
simStruct.params.sigmaBeta   = 0.005; % Estimation; It shall be small so that beta is rarely > 1
simStruct.params.sigmaPolDef = 0.0100; % Arbitrary


% Override with Galí Textbook, Chapter 3
if strcmp(calibrationName, 'textbook')
    
    simStruct.params.rhoA      = 0.9;
    simStruct.params.rhoM      = 0.5;
    simStruct.params.rhoR      = 0.9;
    simStruct.params.rhoD      = 0.5;

    simStruct.params.sigmaA    = 0.01;
    simStruct.params.sigmaM    = 0.0025;
    simStruct.params.sigmaR    = 0.01;
    simStruct.params.sigmaD    = 0.01;

elseif strcmp(calibrationName, 'Brazil')
    
    simStruct.params.rhoA      = 0.91; % Castro (2015)
    simStruct.params.rhoM      = 0.84; % Castro (2015), inflation target shock
    
end

%% MONETARY POLICY RULES
% fixedIntercept_priceLevel,
% rIntercept_priceLevel,
% rRFIntercept_priceLevel,
% optimalIntercept_priceLevel
% fixedIntercept_inflation,
% rIntercept_inflation,
% rRFIntercept_inflation
if endsWith(policyRule, '_inflation')
    simStruct.params.piiBar         = 0.011;    %0.011 Castro (2015)
    simStruct.params.phiNr          = 0; % 0.79 Castro (2015)
    simStruct.params.phi            = 2.43;     % Castro (2015) monetary policy reaction to inflation gap
    simStruct.params.phi_Exp        = 0;        % monetary policy reaction to expected inflation gap
    simStruct.params.phi_Y          = 0;     % 0.16 Castro (2015) monetary policy reaction to output gap
    simStruct.params.phi_dY         = 0;        % monetary policy reaction to output growth
    simStruct.params.elbRate        = 0;        % effective lower bound rate
elseif endsWith(policyRule, '_priceLevel')
    simStruct.params.piiBar         = 0;
    simStruct.params.phiNr          = 0;
    simStruct.params.phi            = 0.5; %monetary policy reaction to price level gap
    simStruct.params.phi_Exp        = 0; %monetary policy reaction to expected price level gap
    simStruct.params.phi_Y          = 0; %monetary policy reaction to output gap
    simStruct.params.phi_dY         = 0; %monetary policy reaction to output growth
    simStruct.params.elbRate        = 0;        % effective lower bound rate
end

% Override with Galí Textbook, Chapter 3
if strcmp(calibrationName, 'textbook')
    
    simStruct.params.piiBar    = 0;
    
    if endsWith(policyRule, '_priceLevel')
        simStruct.params.phi       = 0.5 ;
        simStruct.params.phi_Exp   = 0.0 ;
        simStruct.params.phi_Y     = 0.5 / 4 ;
        simStruct.params.phi_dY    = 0 ;
    elseif endsWith(policyRule, '_inflation')
        simStruct.params.phi       = 1.5 ;
        simStruct.params.phi_Exp   = 0.0 ;
        simStruct.params.phi_Y     = 0.5 / 4 ;
        simStruct.params.phi_dY    = 0 ;
    end
    
    simStruct.params.phiNr    = 0.0; 

elseif strcmp(calibrationName, 'Brazil')
    
    simStruct.params.piiBar    = 0.011; % Castro (2015)
    
    if endsWith(policyRule, '_priceLevel')
        simStruct.params.phiNr     = 0 ;
        simStruct.params.phi       = 0.5;
        simStruct.params.phi_Exp   = 0.0 ;
        simStruct.params.phi_Y     = 0.5 / 4 ;
        simStruct.params.phi_dY    = 0 ;
    elseif endsWith(policyRule, '_inflation')
        simStruct.params.phiNr     = 0.79 ; % Castro (2015)
        simStruct.params.phi       = 0 ;    % Castro (2015)
        simStruct.params.phi_Exp   = 2.43 ; % Castro (2015)
        simStruct.params.phi_Y     = 0.16 ; % Castro (2015)
        simStruct.params.phi_dY    = 0 ;
    end
    
end

%% FISCAL PARAMETERS
if hasGovernment
    simStruct.params.zYBar       = 0.142; % Government Finance Statistics IMF
    simStruct.params.tauBar      = 0.391; % Government Finance Statistics IMF
    simStruct.params.tauMaxBar   = simStruct.params.chi/(1+simStruct.params.chi); % Government Finance Statistics IMF
    simStruct.params.gammaTau    = 0.108; % Estimation with IMF and BCB data
    simStruct.params.deltaBar    = 0.05;  % ??? 20% annual ????, Arbitrary default haircut
    
    if strcmp(govExpenses, 'Actual')
        simStruct.params.gYBar      = 0.206; % Government Finance Statistics IMF
    elseif strcmp(govExpenses, 'High')
        simStruct.params.gYBar      = 0.40; % Arbitrarily high number
    elseif strcmp(govExpenses, 'Zero')
        simStruct.params.gYBar      = 0; % No gov expenses
    end
    
    simStruct.params.debtWedge  = 0;    % Dummy value, it is calculated at the steady state file
else
    simStruct.params.zYBar      = 0;
    simStruct.params.tauBar     = 0;
    simStruct.params.gammaTau   = 0; 
    simStruct.params.deltaBar   = 0;  %default haircut
    
    simStruct.params.gYBar      = 0;
    
    simStruct.params.debtWedge  = 0;
end


if strcmp(calibrationName, 'Brazil')
   
    %%%%%%% GOVERNMENT
    if hasGovernment
        simStruct.params.rhoG           = 0.7947;
        simStruct.params.sigmaG         = 0.0117;
        simStruct.params.gammaTau       = 0.005;
        simStruct.params.tauBar         = 0.3870;
        simStruct.params.gYBar          = 0.1990;
        simStruct.params.zYBar          = 0.1320;
        simStruct.params.bYBar          = 5.0909;
        simStruct.params.deltaBar       = 0;
    end
    %%%%%%%
    
end

%% Find \overline{A} and \eta and tfpLoss
syms aTilde_dummy recTilde_dummy eta_dummy Psi_dummy
eqns = [
        aTilde_dummy      == ( simStruct.params.yBar / ...
            (simStruct.params.kBar*(1-simStruct.params.tauBar) ...
            *(simStruct.params.kBar/eta_dummy)^(1/simStruct.params.chi)) )^(1/(1+1/simStruct.params.chi)) ...
            / ((1-simStruct.params.defBar) + simStruct.params.defBar*simStruct.params.ratioRecA);
        simStruct.params.nBar == ((1-simStruct.params.tauBar)*simStruct.params.kBar/eta_dummy * Psi_dummy )^(1/simStruct.params.chi);
        recTilde_dummy    == simStruct.params.ratioRecA * aTilde_dummy;
        Psi_dummy         == ...
            (1-simStruct.params.defBar)*aTilde_dummy ...
            + simStruct.params.defBar*recTilde_dummy;
        ];
soln = vpasolve(eqns, symvar(eqns));
simStruct.params.aTildeBar      = double(soln.aTilde_dummy);
simStruct.params.recTildeBar    = double(soln.recTilde_dummy);
simStruct.params.eta            = double(soln.eta_dummy);
simStruct.params.PsiBar         = double(soln.Psi_dummy);
clear('aTilde_dummy');
clear('recTilde_dummy');
clear('eta_dummy');
clear('Psi_dummy');

% TFP Loss
% https://ftp.ibge.gov.br/Contas_Nacionais/Sistema_de_Contas_Nacionais/2018/tabelas_xls/PIBReal_1947_2018.xls
simStruct.params.tfpLoss     = 0.0238;  % 0.0238 TFP shock is consistent with -4.3% YoY, which is the GDP growth of Brazil in 1990

%% FISCAL PARAMETERS 2

if hasGovernment
   
    %%%%%%%%%%%%%%% DEFAULT
    simStruct.params.fiscalLimit_mu     = fiscalLimit_mu;
    simStruct.params.fiscalLimit_std    = fiscalLimit_std;
    
    if govAccumDebt == false
        simStruct.params.bYBar      = 0; % No debt
        
        % Parameterize the reduced form version of the fiscal limits
        approximateFiscalLimits;
    elseif strcmp(debtLevel, 'Actual')
        simStruct.params.bYBar      = 0.619*4; % General Government Gross Debt - BCB
        
        % Parameterize the reduced form version of the fiscal limits
        approximateFiscalLimits;
    elseif strcmp(debtLevel, 'Zero')
        simStruct.params.bYBar      = 0; % No debt
        
        % Parameterize the reduced form version of the fiscal limits
         approximateFiscalLimits;
    elseif strcmp(debtLevel, 'High')
        
        if ~exist('defPolTarget', 'var')
            defPolTarget = 0.0025; % Target of quarterly default probability, [] for to be calibrated
        end
        simStruct.params.bYBar = 0.619*4; % Start candidate = General Government Gross Debt - BCB

        % Parameterize the reduced form version of the fiscal limits
        approximateFiscalLimits;

        calibrateDebtToDefaultProb;
        simStruct.params.bYBar = debtTargetDefProb;
        %simStruct.params.gammaTau = gammaTauTargetDefProb;
        
        % Again!!
        approximateFiscalLimits;


        %simStruct.params.bYBar = fzero(@(bYTarget) defPolTarget - fFiscalLimit_defProb(...
        %    bYTarget * simStruct.params.yBar , ...
        %    0, ...
        %    0, ...
        %    0, ...
        %    0, ...
        %    0), 1);
        
    end

%     defProb = fFiscalLimit_defProb(...
%                 simStruct.params.bYBar * simStruct.params.yBar , ...
%                 simStruct.params.aTildeBar, ...
%                 simStruct.params.recTildeBar, ...
%                 simStruct.params.defBar, ...
%                 simStruct.params.gYBar * simStruct.params.yBar, ...
%                 simStruct.params.beta);
    defProb = fFiscalLimit_defProb(...
                simStruct.params.bYBar * simStruct.params.yBar , ...
                0, ...
                0, ...
                0, ...
                0, ...
                0);
    deltaExp    = defProb * simStruct.params.deltaBar;
    rPolicy     = -1 + 1/(simStruct.params.beta * (1 - deltaExp));

    
    %%%%%%%%%%%%%% CALCULATE INTEREST RATE TO CALCULATE LUMP SUM TAXES
    rRN         = -1 + 1/(simStruct.params.beta);
    nr          = (1 + rRN) * (1+simStruct.params.piiBar) - 1;

    if startsWith(policyRule, "fixedIntercept_") || startsWith(policyRule, "rRN_") || startsWith(policyRule, "rRF_")  || startsWith(policyRule, "rEF_")
        iota        = rRN ;
        nrPolicy          = (1 + iota) * (1+simStruct.params.piiBar) - 1;
    elseif  startsWith(policyRule, "rRNwithRiskFree_")
        %iota        = rPolicy ;
        nrPolicy          = (1 + rRN) * (1+simStruct.params.piiBar) - 1;
    elseif startsWith(policyRule, "polDefAdjusted_")
        iota = rPolicy;
        nrPolicy = (1 + iota) * (1+simStruct.params.piiBar) - 1;
    end

    if endsWith(policyRule, "_priceLevel")
        pLevel               = simStruct.params.pLevelBar;
    end
    
     %%%%%%%%%%% LUMP-SUM TAX CALIBRATION
    syms debtWedge_sym
    deltaBar_NoDefault = 0;
    eqns = [
                simStruct.params.bYBar == ...
                (1 + nrPolicy) * (simStruct.params.gYBar + simStruct.params.zYBar + debtWedge_sym - simStruct.params.tauBar) ...
                / (1  -(1 + nrPolicy)*(1-deltaBar_NoDefault)/(1+simStruct.params.piiBar));
            ];
    soln            = vpasolve(eqns, symvar(eqns));
    simStruct.params.debtWedge  = double(soln);
    %%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    simStruct.params.fiscalLimit_mu     = 0;
    simStruct.params.fiscalLimit_std    = 0;
    
    simStruct.params.bYBar      = 0;
    
end


if hasGovernment
    
    if strcmp(calibrationName, 'Brazil')
        simStruct.params.bYBar          = 5.0909;
    end
    
end
    
%% RISE PARAMETERS
simStruct.params.ytrend             = 0;
simStruct.params.polDefWedge        = 0;     % difference at the steady-state between the default prob. in FX and LC
simStruct.params.unempSS            = simStruct.params.nBar;
simStruct.params.sigmaK             = 0.01;

simStruct.params.stderr_obs_y                       = 0.01;
simStruct.params.stderr_obs_c                       = 0.01;
simStruct.params.stderr_obs_g                       = 0.01;
simStruct.params.stderr_obs_unempRate               = 0.01;
simStruct.params.stderr_obs_wage                    = 0.01;
simStruct.params.stderr_obs_pii                     = 0.003;
simStruct.params.stderr_obs_defLC                   = 0.0025;
simStruct.params.stderr_obs_nr                      = 0.0025;
simStruct.params.stderr_obs_swap_PreDI_3m           = 0.0025;
simStruct.params.stderr_obs_swap_PreDI_6m           = 0.0025;
simStruct.params.stderr_obs_pii_FOCUS_Median        = 0.005;
simStruct.params.stderr_obs_y_FOCUS_Median_BOQ      = 0.01;
simStruct.params.stderr_obs_y_FOCUS_0               = 0.01;
simStruct.params.stderr_obs_y_FOCUS_25              = 0.01;
simStruct.params.stderr_obs_y_FOCUS_75              = 0.01;
simStruct.params.stderr_obs_y_FOCUS_100             = 0.01;
simStruct.params.stderr_obs_y_FOCUS_20_EOQ          = 0.01;
simStruct.params.stderr_obs_y_FOCUS_80_EOQ          = 0.01;

%% Load estimation results

% Parameters that enter into the calculation of fiscal limits
%aTildeSS, recTildeSS, gSS, defSS, ySS, zSS, tLSSS, yMaxSS, tauMaxSS, ...
%beta, chi, eta, kkSS, sigma, alphaG, gammaGPSI, fracNR, ...
%rhoA, rhoR, rhoD, rhoGG, rhoGY, ...
%sigmaA, sigmaR, sigmaD, sigmaG, ...

% NOTE THAT ONLY DEEP PARAMETERS ARE CHANGED
% DEPENDENT PARAMETERS HAVE TO BE UPDATED MANUALLY BELOW
if ~strcmp(loadEstimParameters, 'None')
    load('estimPosteriorResults');
    
    paramEstimNames = fieldnames(estimPostSimData);
    
    if strcmp(loadEstimParameters, 'OnlyShocks')
        paramEstimNames = { 'rhoA'; 'rhoBeta'; 'rhoGG'; 'rhoGY'; 'rhoPolDef'; ...
                            'sigmaA'; 'sigmaBeta'; 'sigmaG'; 'sigmaM'; 'sigmaPolDef'};
    end
    
    % Posterior mean
    for iParam = 1:length(paramEstimNames)
        simStruct.params.(paramEstimNames{iParam}) = estimPostSimData.(paramEstimNames{iParam}).mean_sim;
    end
    
end

%% Turn off shocks

for iShock = 1:length(shocksTurnedOff)
    simStruct.params.(shocksTurnedOff(iShock)) = 0;
end

%% Override estimation results

% Check whether beta standard deviation is too large
if simStruct.params.beta + 3*simStruct.params.sigmaBeta > 1
    %simStruct.params.beta        = 0.989;
    %simStruct.params.sigmaBeta  = 0.001;
    %simStruct.params.sigma      = 1.3;
    %simStruct.params.gammaGPSI  = 0.2;
    %simStruct.params.chi        = 1.3;
end

% Check whether policy rule is simplified
if simplifyPolicyRule
    simStruct.params.piiBar         = 0.011;
    simStruct.params.phiNr          = 0;
    simStruct.params.phi_Exp        = 0; %monetary policy reaction to expected price level gap
    simStruct.params.phi_Y          = 0; %monetary policy reaction to output gap
    simStruct.params.phi_dY         = 0; %monetary policy reaction to output growth
end