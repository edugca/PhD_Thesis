%% Steady state values for REGIME 1

bindDef = 0;
bindTax = 0;
bindFis = 0;

%% BASIC VARIABLES

simStruct.ss.pLevelSS   = simStruct.params.pLevelBar;
simStruct.ss.piiSS      = simStruct.params.piiBar;
simStruct.ss.piiExpSS   = simStruct.params.piiBar;

simStruct.ss.aTildeSS   = simStruct.params.aTildeBar;
simStruct.ss.recTildeSS = simStruct.params.recTildeBar;

if strcmp(defRProcessType, "Null")
    simStruct.ss.defSS     = 0 ;
    simStruct.ss.defExpSS  = 0 ;
elseif strcmp(defRProcessType, "MS_Exogenous") || strcmp(defRProcessType, "MS_Endogenous")
    simStruct.ss.defSS     = bindDef ;
    simStruct.ss.defExpSS  = bindDef ;
else
    simStruct.ss.defSS     = simStruct.params.defBar ;
    simStruct.ss.defExpSS  = simStruct.params.defBar ;
end

simStruct.ss.PsiSS      = (1. - simStruct.ss.defSS) * simStruct.ss.aTildeSS ...
                            + simStruct.ss.defSS * simStruct.ss.recTildeSS;
simStruct.ss.PsiExpSS   = simStruct.ss.PsiSS;

simStruct.ss.kkSS         = simStruct.params.kBar;

simStruct.ss.mShockSS       = 1;
simStruct.ss.mShockExpSS    = simStruct.ss.mShockSS;

simStruct.ss.bbetaSS    = simStruct.params.beta;

%% FISCAL VARIABLES (PART 1)

% Check whether the economy has government
if hasGovernment == false
    simStruct.ss.zYSS   = 0;
    simStruct.ss.gYSS   = 0;
    simStruct.ss.bYSS   = 0;
    simStruct.ss.tauSS  = 0;
else
    simStruct.ss.zYSS       = simStruct.params.zYBar;
    simStruct.ss.gYSS       = simStruct.params.gYBar;
    simStruct.ss.bYSS       = simStruct.params.bYBar;
    simStruct.ss.tauSS      = simStruct.params.tauBar;
end

%% REAL SECTOR (PART 1)

simStruct.ss.ySS        = simStruct.params.yBar;
simStruct.ss.nSS        = simStruct.params.nBar;
simStruct.ss.wpSS       = kBar .* simStruct.ss.PsiExpSS;

%% FISCAL VARIABLES (PART 2)

simStruct.ss.taxSS      = simStruct.ss.tauSS * simStruct.ss.ySS;

simStruct.ss.gSS        = simStruct.ss.gYSS * simStruct.ss.ySS;
simStruct.ss.zSS        = simStruct.ss.zYSS * simStruct.ss.ySS;

simStruct.ss.deltaSS    = 0; % NOT DEFAULT at the steady state
simStruct.ss.bSS        = simStruct.ss.bYSS * simStruct.ss.ySS;
simStruct.ss.bdSS       = (1 - simStruct.ss.deltaSS) * simStruct.ss.bSS;

if govBondsAreRiskFree == false
    
    simStruct.ss.muFisLimSS = fFiscalLimit_mu(...
                                simStruct.ss.aTildeSS, ...
                                simStruct.ss.recTildeSS, ...
                                simStruct.ss.defSS, ...
                                simStruct.ss.gSS, ...
                                simStruct.ss.bbetaSS);
    simStruct.ss.stdFisLimSS  = fFiscalLimit_std(...
                                simStruct.ss.aTildeSS, ...
                                simStruct.ss.recTildeSS, ...
                                simStruct.ss.defSS, ...
                                simStruct.ss.gSS, ...
                                simStruct.ss.bbetaSS);
                            
    simStruct.ss.probDefFisLimSS  = fFiscalLimit_defProb(...
                                simStruct.ss.bSS , ...
                                simStruct.ss.aTildeSS, ...
                                simStruct.ss.recTildeSS, ...
                                simStruct.ss.defSS, ...
                                simStruct.ss.gSS, ...
                                simStruct.ss.bbetaSS);

    simStruct.ss.polDefSS   = fFiscalLimit_defProb(...
                                simStruct.ss.bSS , ...
                                simStruct.ss.aTildeSS, ...
                                simStruct.ss.recTildeSS, ...
                                simStruct.ss.defSS, ...
                                simStruct.ss.gSS, ...
                                simStruct.ss.bbetaSS);
                            
else
   
   simStruct.ss.muFisLimSS  = 0; 
   simStruct.ss.stdFisLimSS  = 0;
   simStruct.ss.probDefFisLimSS  = 0; 
   simStruct.ss.polDefSS    = 0;
   
end
                        
% %%%% NOT BEING USED
% fiscalLimit = 0;
% for j=1:nFiscalLimit
%     fiscalLimit = fiscalLimit + (beta.^(j-1) .* (simStruct.ss.taxMaxSS - simStruct.ss.gSS - simStruct.ss.zSS));
% end
% % Check whether there is fiscal limit data
% if sum(isnan(vec(fiscalLimit_mu))) == length(vec(fiscalLimit_mu))
%     % No data
%     disp('There is no fiscal limit data!');
%     pdf_delta = makedist('Normal');
% else
%     % There is data
%     pdf_delta = pickFiscalLimit(simStruct, fiscalLimit_mu, fiscalLimit_std, 0,0,0,0);
% end
% f_delta = @(x) cdf(pdf_delta,x);
% simStruct.ss.polDefSS   = f_delta(simStruct.ss.bSS);

%% REAL SECTOR (PART 2)

simStruct.ss.cSS        = simStruct.ss.ySS - simStruct.ss.gSS;
simStruct.ss.cNRSS      = (1-simStruct.ss.tauSS) * simStruct.ss.wpSS * simStruct.ss.nSS ;
simStruct.ss.cRSS       = 1/(1-simStruct.params.fracNR) * (simStruct.ss.cSS - simStruct.params.fracNR * simStruct.ss.cNRSS) ; 
simStruct.ss.uCSS       = ( simStruct.ss.cRSS + simStruct.params.alphaG*simStruct.ss.gSS - simStruct.params.eta.*simStruct.ss.nSS.^(1+simStruct.params.chi)./(1+simStruct.params.chi) ).^(-simStruct.params.sigma);

simStruct.ss.cStarSS        = simStruct.ss.ySS - (1-simStruct.params.alphaG)*simStruct.ss.gSS;
simStruct.ss.cNRStarSS      = (1-simStruct.ss.tauSS) * simStruct.ss.wpSS * simStruct.ss.nSS + simStruct.params.alphaG*simStruct.ss.gSS ;
simStruct.ss.cRStarSS       = 1/(1-simStruct.params.fracNR) * (simStruct.ss.cSS - simStruct.params.fracNR * simStruct.ss.cNRSS) ; 

simStruct.ss.welfareR    = 1/(1-simStruct.ss.bbetaSS) * 1/(1-simStruct.params.sigma)*(simStruct.ss.cRSS - simStruct.params.eta*simStruct.ss.nSS^(1+simStruct.params.chi)/(1+simStruct.params.chi))^(1-simStruct.params.sigma);
simStruct.ss.welfareNR   = 1/(1-simStruct.ss.bbetaSS) * 1/(1-simStruct.params.sigma)*(simStruct.ss.cNRSS - simStruct.params.eta*simStruct.ss.nSS^(1+simStruct.params.chi)/(1+simStruct.params.chi))^(1-simStruct.params.sigma);

%% MONETARY POLICY and PORTFOLIO

simStruct.ss.rRNSS        = 1/simStruct.params.beta - 1;
simStruct.ss.rPolicySS    = 1/(simStruct.params.beta * (1-simStruct.ss.deltaSS)) - 1;
%simStruct.ss.rGapSS     = simStruct.ss.rPolicySS - simStruct.ss.iotaSS ;
simStruct.ss.nrSS         = (1 + simStruct.ss.rPolicySS) .* (1+simStruct.ss.piiSS) - 1;

%% LUMP-SUM TAX CALIBRATION

syms debtWedge_sym
eqns = [
            simStruct.ss.bSS == (1 + simStruct.ss.nrSS) * (simStruct.ss.gSS + simStruct.ss.zSS + debtWedge_sym - simStruct.ss.taxSS) / (1  -(1 + simStruct.ss.nrSS)*(1-simStruct.ss.deltaSS)/(1+simStruct.ss.piiSS));
        ];
soln            = vpasolve(eqns, symvar(eqns));
simStruct.ss.tLSSS  = double(soln);
clear('debtWedge_sym');

%% AUXILIARY EXPECTATION VARIABLES

simStruct.ss.aTildeExpSS    = simStruct.ss.aTildeSS;
simStruct.ss.recTildeExpSS  = simStruct.ss.recTildeSS;
simStruct.ss.defExpSS       = simStruct.ss.defSS;
simStruct.ss.gExpSS         = simStruct.ss.gSS;
simStruct.ss.mShockExpSS    = simStruct.ss.mShockSS;
simStruct.ss.zExpSS         = simStruct.ss.zSS;
simStruct.ss.PsiExpSS       = simStruct.ss.PsiSS;

simStruct.ss.tauExpSS       = simStruct.ss.tauSS;
simStruct.ss.taxExpSS       = simStruct.ss.taxSS;
simStruct.ss.yExpSS         = simStruct.ss.ySS;
simStruct.ss.cExpSS         = simStruct.ss.cSS;
simStruct.ss.cRExpSS        = simStruct.ss.cRSS;
simStruct.ss.cNRExpSS       = simStruct.ss.cNRSS;
simStruct.ss.nExpSS         = simStruct.ss.nSS;
simStruct.ss.uCExpSS        = simStruct.ss.uCSS;

%% MAXIMUM VARIABLES

if hasGovernment

    simStruct.ss.tauMaxSS          = simStruct.params.chi/(1+simStruct.params.chi);
    simStruct.ss.taxMaxSS          = simStruct.ss.tauMaxSS * simStruct.ss.kkSS * ((1-simStruct.ss.tauMaxSS) * simStruct.ss.kkSS/simStruct.params.eta * simStruct.ss.PsiExpSS)^(1/simStruct.params.chi) * simStruct.ss.PsiSS;
    simStruct.ss.wpMaxSS           = simStruct.ss.kkSS * simStruct.ss.PsiExpSS;
    simStruct.ss.nMaxSS            = ( (1-simStruct.ss.tauMaxSS) * simStruct.ss.kkSS/simStruct.params.eta * simStruct.ss.PsiExpSS )^(1/simStruct.params.chi);
    simStruct.ss.yMaxSS            = simStruct.ss.kkSS * ( (1-simStruct.ss.tauMaxSS) * simStruct.ss.kkSS/simStruct.params.eta * simStruct.ss.PsiExpSS )^(1/simStruct.params.chi) * simStruct.ss.PsiSS;
    simStruct.ss.cMaxSS            = simStruct.ss.yMaxSS - simStruct.ss.gSS;
    simStruct.ss.cNRMaxSS          = (1-simStruct.ss.tauMaxSS) * simStruct.ss.wpMaxSS * simStruct.ss.nMaxSS ;
    simStruct.ss.cRMaxSS           = 1/(1-simStruct.params.fracNR) * (simStruct.ss.cMaxSS - simStruct.params.fracNR * simStruct.ss.cNRMaxSS) ; 
    simStruct.ss.uCMaxSS           = ( simStruct.ss.cRMaxSS + simStruct.params.alphaG*simStruct.ss.gSS - simStruct.params.eta*simStruct.ss.nMaxSS^(1+simStruct.params.chi)/(1+simStruct.params.chi) )^(-simStruct.params.sigma) ;
    
else
    
    simStruct.ss.tauMaxSS          = simStruct.ss.tauSS;
    simStruct.ss.taxMaxSS          = simStruct.ss.taxSS;
    simStruct.ss.wpMaxSS           = simStruct.ss.wpSS;
    simStruct.ss.nMaxSS            = simStruct.ss.nSS;
    simStruct.ss.yMaxSS            = simStruct.ss.ySS;
    simStruct.ss.cMaxSS            = simStruct.ss.cSS;
    simStruct.ss.cRMaxSS           = simStruct.ss.cRSS;
    simStruct.ss.cNRMaxSS          = simStruct.ss.cNRSS;
    simStruct.ss.uCMaxSS           = simStruct.ss.uCSS;
    
end