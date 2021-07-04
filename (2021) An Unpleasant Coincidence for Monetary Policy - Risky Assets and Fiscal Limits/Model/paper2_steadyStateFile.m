function [ss,newp,retcode] = paper2_steadyStateFile(m, ss, p, d, id)

% Args:
%    m (rise | dsge): model object (not always needed)
%    ss (vector): endo_nbr x 1 vector of initial steady state
%    p (struct): parameter structure
%    d (struct): definitions
%    id (vector): location of the variables to calculate

% Returns:
%    :
%
%    - **ss** []: endo_nbr x 1 vector of updated steady state
%    - **newp** [struct]: structure containing updated parameters if any
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0

%% PRELIMINARY
retcode = 0;
if nargin == 1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    ss = m.endogenous.name;

    % list of parameters to be computed during steady state calculation
    %-------------------------------------------------------------------
    newp={'ySS', 'debtWedge'};
    
    return
    
end

% no parameters to update or create in the steady state file
%-----------------------------------------------------------
newp=[];

%% BASIC VARIABLES

aTilde          = p.aTildeBar;
aTildeExp       = p.aTildeBar;
recTilde        = p.recTildeSS;
recTildeExp     = recTilde;

Psi             = aTilde ;
PsiExp          = aTildeExp ;

kk              = p.kBar;

mShock          = p.mShockSS;
mShockExp       = mShock;

bbeta           = p.beta;

%% FISCAL VARIABLES (PART 1)

if m.user_data.conf_hasGovernment
    
    tau             = (1-p.bindTax)*p.tauBar + p.bindTax*p.tauMaxBar;
    
    b  = p.bYBar * p.yBar ;
    
    if strcmp(m.user_data.conf_approximateDefaultProb, 'logistic')
        probDefFisLim   = 1/(1 ...
                                + exp( p.probDefFisLim_Param_0  ...
                                + p.probDefFisLim_Param_B*(b - p.muFisLim_Intercept) ));
    elseif strcmp(m.user_data.conf_approximateDefaultProb, 'linear')
        probDefFisLim   = p.probDefFisLim_Param_0 + p.probDefFisLim_Param_B*(b - p.muFisLim_Intercept) ;
    end
    
    if m.user_data.conf_govBondsAreRiskFree
        % Set fiscal limits to infinity
        ddelta = 0 ;
    else
        % Default?
        ddelta = p.bindFis * p.deltaBar ;
        %ddelta = deltaBar ;
    end
    
else
    
    tau     = 0;
    b       = 0;
    probDefFisLim = 0;
    ddelta   = 0;
    
end

ddeltaExp    = ddelta;

%% (FLEXIBLE/STICKY) PRICES / REAL SECTOR (PART 1)

mc          = (p.thetaElast-1)/p.thetaElast /(1-p.labSub);
adjCosts    = 0;

n               = ( (mc*aTilde*kk)/(p.eta/(1-tau)) )^(1/p.chi);
%n               = ( wp * (1-tau) / p.eta )^(1/p.chi);
wp              = (1/(1-tau)) * p.eta * (n^p.chi) ;
y               = aTilde * kk * n;

tax             = tau*(n*wp) + tau*(y - n*wp - adjCosts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n               = ( (1-tau) * kk/p.eta * aTilde )^(1/p.chi);
% wp              = (1/(1-tau)) * p.eta * n^p.chi ;
% y               = kk * ( (1-tau) * kk/p.eta * aTilde )^(1/p.chi) * aTilde;
% 
% tax             = tau * kk * ((1-tau) * kk/p.eta * aTilde)^(1/p.chi) * aTilde;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LAG_1_y         = y;
newp.ySS        = y;

%% MONETARY POLICY

rRN             = -1 + 1/(bbeta);
rNa             = -1 + 1/(bbeta);
rEF             = -1 + 1/(bbeta);

%rPolicy         = -1 + 1 / ( (1 - p.bindFis) * bbeta +  ( p.bindFis*(1-p.deltaBar) ) * bbeta ) ;
rPolicy         = -1 + 1 / ( (1 - probDefFisLim) * bbeta +  ( probDefFisLim*(1-p.deltaBar) ) * bbeta ) ;

Pii             = p.bindELB * ( (1+p.elbRate)/(1+rPolicy) ) + (1 - p.bindELB) * (1 + p.piiBar);
LAG_1_Pii       = Pii;
LAG_2_Pii       = Pii;
LAG_3_Pii       = Pii;

Pii_YoY         = Pii^4;
LEAD_1_Pii_YoY  = Pii_YoY;
LEAD_2_Pii_YoY  = Pii_YoY;
LEAD_3_Pii_YoY  = Pii_YoY;

if startsWith(m.user_data.conf_policyRule, "fixedIntercept_") || startsWith(m.user_data.conf_policyRule, "rRN_") || startsWith(m.user_data.conf_policyRule, "polDefAdjusted_")
    iota = rPolicy;
    nrPolicy    = - 1 + p.bindELB*(1+p.elbRate) + (1-p.bindELB) * ( (1 + iota) * Pii );
    nrPolicyExp = - 1 + (1 + nrPolicy)^p.phiNr * ( (1 + iota) * (1 + p.piiBar) * (Pii / (1 + p.piiBar))^p.phi )^(1-p.phiNr);
elseif startsWith(m.user_data.conf_policyRule, "rRF_") || startsWith(m.user_data.conf_policyRule, "rRNwithRiskFree_") || startsWith(m.user_data.conf_policyRule, "rEF_") || startsWith(m.user_data.conf_policyRule, "rNa_")
    iota        = rRN ;
    nrPolicy    = - 1 + p.bindELB*(1+p.elbRate) + (1-p.bindELB) * ( (1 + iota) * Pii );
    nrPolicyExp = - 1 + (1 + nrPolicy)^p.phiNr * ( (1 + iota) * (1 + p.piiBar) * (Pii / (1 + p.piiBar))^p.phi )^(1-p.phiNr);
end

if endsWith(m.user_data.conf_policyRule, "_priceLevel")
    pLevel               = p.pLevelBar;
end

%% FISCAL VARIABLES (PART 2)

if m.user_data.conf_hasGovernment
    
    g  = p.gYBar * p.yBar;
    z  = p.zYBar * p.yBar;
    gY = g/y;
    bY = b/y;
    %newp.bYSS = bY;
    
    if m.user_data.conf_govBondsAreRiskFree
        muFisLim        = 0;
        stdFisLim       = 0;
        probDefFisLim   = 0;
        drawnFisLim     = 0;
        
        polDef      = 0;
        p.polDefSS  = polDef;
        rPolicy     = -1 + 1/p.beta;
        rGov        = rPolicy;
        nrGov       = - 1 + p.bindELB*(1 + p.elbRate) + (1 - p.bindELB) * ( (1 + rGov) * Pii ) ;
    else
        % Steady state values of muFisLim and stdFisLim
        %rise_fiscalLimit;
        muFisLim        = p.muFisLim_Intercept;
        stdFisLim       = p.stdFisLim_Intercept;
        drawnFisLim     = muFisLim;

        polDef      = p.bindFis;
        p.polDefSS  = polDef;
        %rGov        = -1 + 1 / ( (1 - polDef) * p.beta +  ( polDef*(1-p.deltaBar) ) * p.beta ) ;
        rGov        = -1 + 1 / ( (1 - probDefFisLim) * p.beta +  ( probDefFisLim*(1-p.deltaBar) ) * p.beta ) ;
        nrGov       = - 1 + p.bindELB*(1 + p.elbRate) + (1 - p.bindELB) * ( (1 + rGov) * Pii ) ;
        
%         muFisLim = fFiscalLimit_mu(...
%                                 aTilde, ...
%                                 recTilde, ...
%                                 def, ...
%                                 g);
%                             
%         stdFisLim = fFiscalLimit_std(...
%                                 aTilde, ...
%                                 recTilde, ...
%                                 def, ...
%                                 g);
    end
    
    if startsWith(m.user_data.conf_policyRule, "polDefAdjusted_")
        %rPolicy            = -1 + 1 / ( (1 - polDef) * p.beta +  ( polDef*(1-p.deltaBar) ) * p.beta ) ;
        rPolicy             = rGov ;
        iota                = rPolicy;
        nrPolicy            = - 1 + p.bindELB*(1+p.elbRate) + (1-p.bindELB) * ( (1 + iota) * Pii ) ;
        nrPolicyExp         = - 1 + (1 + nrPolicy)^p.phiNr * ( (1 + iota) * (1 + p.piiBar) * (Pii / (1 + p.piiBar))^p.phi )^(1-p.phiNr) ;
    elseif startsWith(m.user_data.conf_policyRule, "fixedIntercept_") || startsWith(m.user_data.conf_policyRule, "rRN_") || startsWith(m.user_data.conf_policyRule, "rEF_") || startsWith(m.user_data.conf_policyRule, "rNa_")
        %rPolicy             = -1 + 1 / ( (1 - polDef) * p.beta +  ( polDef*(1-p.deltaBar) ) * p.beta ) ;
        rPolicy             = rGov ;
        iota                = rPolicy ;
        %nrPolicy            = (1 + iota) * Pii - 1;
        nrPolicy            = - 1 + p.bindELB*(1+p.elbRate) + (1-p.bindELB) * ( (1 + iota) * Pii );
        nrPolicyExp         = - 1 + (1 + nrPolicy)^p.phiNr * ( (1 + iota) * (1 + p.piiBar) * (Pii / (1 + p.piiBar))^p.phi )^(1-p.phiNr) ;
    elseif startsWith(m.user_data.conf_policyRule, "rRNwithRiskFree_")
        %rPolicy             = -1 + 1 / ( (1 - polDef) * p.beta +  ( polDef*(1-p.deltaBar) ) * p.beta ) ;
        rPolicy             = -1 + 1 / p.beta ;
        iota                = rRN ;
        %nrPolicy            = (1 + rRN) * Pii - 1;
        nrPolicy            = - 1 + p.bindELB*(1+p.elbRate) + (1-p.bindELB) * ( (1 + rRN) * Pii );
        nrPolicyExp         = - 1 + (1 + nrPolicy)^p.phiNr * ( (1 + rRN) * (1 + p.piiBar) * (Pii / (1 + p.piiBar))^p.phi )^(1-p.phiNr) ;
    elseif startsWith(m.user_data.conf_policyRule, "rRF_")
        rPolicy             = -1 + 1 / p.beta ;
        iota                = rRN ;
        %nrPolicy            = (1 + rRN) * Pii - 1;
        nrPolicy            = - 1 + p.bindELB*(1+p.elbRate) + (1-p.bindELB) * ( (1 + rRN) * Pii );
        nrPolicyExp         = - 1 + (1 + nrPolicy)^p.phiNr * ( (1 + rRN) * (1 + p.piiBar) * (Pii / (1 + p.piiBar))^p.phi )^(1-p.phiNr) ;
    end
    
    %%%%%%%%%%% LUMP-SUM TAX CALIBRATION
    syms debtWedge_sym
    eqns = [
                b == (1 + nrGov) * (g + z + debtWedge_sym - tax) / (1  -(1 + nrGov)*(1-ddelta)/Pii);
            ];
    soln            = vpasolve(eqns, symvar(eqns));
    newp.debtWedge  = double(soln);
    %newp.debtWedge  = -0.041948220706830859649508400977328;
    %%%%%%%%%%%
    
    bd              = b ;
    bLag1           = b ;
    
else
    
    g       = 0;
    gY      = 0;
    z       = 0;
    
    bd      = 0;
    bLag1   = 0;
    bY      = 0;
    
    probDefFisLim   = 0;
    muFisLim        = 0;
    stdFisLim       = 0;
    drawnFisLim     = muFisLim;
    polDef          = 0;
    
    rPolicy     = -1 + 1/p.beta;
end

rGap        = rPolicy - rRN ;
rGapIota    = iota - rRN ;


%% REAL SECTOR (PART 2)

c               = y - g;
cNR             = (1-tau) * wp * n + z + newp.debtWedge ;
cR              = 1/(1-p.fracNR) * (c - p.fracNR * cNR) ;

% Natural
nNa     = ( (p.thetaElast - 1)/p.thetaElast * (1-tau)/p.eta * kk * aTilde  )^(1/p.chi) ;
cRNa    = 1/(1-p.fracNR) * ( ( (p.thetaElast - 1)/p.thetaElast * (1-tau)/p.eta )^(1/p.chi) * (kk*aTilde)^(1+1/p.chi) - g - p.fracNR*(( (p.thetaElast - 1)/p.thetaElast * 1/p.eta)^(1/p.chi) * ((1-tau)*kk*aTilde)^(1+1/p.chi) + z ) ) ;

% Efficient
nEF     = ( (1-tau)/p.eta * kk * aTilde  )^(1/p.chi) ;
cREF    = 1/(1-p.fracNR) * ( ( (1-tau)/p.eta )^(1/p.chi) * (kk*aTilde)^(1+1/p.chi) - g - p.fracNR*((1/p.eta)^(1/p.chi) * ((1-tau)*kk*aTilde)^(1+1/p.chi) + z ) ) ;


cStar           = y - (1-p.alphaG)*g ; 
cNRStar         = (1-tau) * wp * n + z + newp.debtWedge + p.alphaG*g;
cRStar          = 1/(1-p.fracNR) * (cStar - p.fracNR * cNRStar) ;

welfareR    = 1/(1-bbeta) * 1/(1-p.sigma)*(cR + p.alphaG*g - p.eta*n^(1+p.chi)/(1+p.chi))^(1-p.sigma);
welfareNR   = 1/(1-bbeta) * 1/(1-p.sigma)*(cNR + p.alphaG*g - p.eta*n^(1+p.chi)/(1+p.chi))^(1-p.sigma);
welfare     = (1-p.fracNR)*welfareR + p.fracNR*welfareNR;

%uC              = ( c - p.eta*n^(1+p.chi)/(1+p.chi) )^(-p.sigma);
uC              = ( cR + p.alphaG*g - p.eta*n^(1+p.chi)/(1+p.chi) )^(-p.sigma);

%yGood           = kk * ( (1-tau) * kk/p.eta * Psi )^(1/p.chi) * aTilde;
%yBad            = kk * ( (1-tau) * kk/p.eta * Psi )^(1/p.chi) * recTilde;

%yGoodExp        = kk * ( (1-tau) * kk/p.eta * Psi )^(1/p.chi) * aTildeExp;
%yBadExp         = kk * ( (1-tau) * kk/p.eta * Psi )^(1/p.chi) * recTildeExp;

%% MAXIMUM VARIABLES

if m.user_data.conf_hasGovernment

    tauMax          = p.tauMaxBar;
    taxMax          = tauMax * kk * ((1-tauMax) * kk/p.eta * Psi)^(1/p.chi) * Psi;
    wpMax           = kk * Psi;
    nMax            = ( (1-tauMax) * kk/p.eta * Psi )^(1/p.chi);
    yMax            = kk * ( (1-tauMax) * kk/p.eta * Psi )^(1/p.chi) * Psi;
    cMax            = yMax - g;
    cNRMax          = (1-tauMax) * wpMax * nMax + z + newp.debtWedge;
    cRMax           = 1/(1-p.fracNR) * (cMax - p.fracNR * cNRMax) ; 
    %uCMax           = ( cMax - p.eta*nMax^(1+p.chi)/(1+p.chi) )^(-p.sigma) ;
    uCMax           = ( cRMax + p.alphaG*g - p.eta*nMax^(1+p.chi)/(1+p.chi) )^(-p.sigma) ;
    
else
    
    tauMax          = tau;
    taxMax          = tax;
    wpMax           = wp;
    nMax            = n;
    yMax            = y;
    cMax            = c;
    cRMax           = cR;
    cNRMax          = cNR;
    uCMax           = uC;
    
end

%% PORTFOLIO VARIABLES

%rGood       = 1/bbeta * ( aTildeExp / Psi ) - 1;
%rGood        = aTildeExp / (p.beta * Psi) - 1;
%nrGood       = 1/bbeta * ( aTildeExp / Psi ) * Pii - 1;

%rBad        = 1/bbeta * ( recTildeExp / Psi ) - 1;
%rBad         = recTildeExp / (p.beta * Psi) - 1;
%nrBad        = 1/bbeta * ( recTildeExp / Psi ) * Pii - 1;

%rGoodOnly   = 1/bbeta * ( aTildeExp / aTilde ) - 1;
%rGoodOnly    = aTildeExp / (p.beta * aTilde) - 1;
%nrGoodOnly   = 1/bbeta * ( aTildeExp / aTilde ) * Pii - 1;

%% AUXILIARY EXPECTATION VARIABLES

tauSdw          = p.tauBar; % it is the shadow of tau

yExp            = y;
cExp            = c;
cRExp           = cR;
cNRExp          = cNR;
nExp            = n;
uCExp           = uC;
gExp            = g;
zExp            = z;
taxExp          = p.tauBar;
taxExp          = tax;

%% ESTIMATION VARIABLES

obs_y                   = p.ytrend;
obs_c                   = p.ytrend;
obs_g                   = p.ytrend;
obs_defLC               = p.polDefWedge;
obs_pii                 = log(Pii);
obs_pii_FOCUS_Median    = log(Pii);
obs_nr                  = nrPolicy;
obs_y_FOCUS_Median_BOQ  = p.ytrend;
obs_y_FOCUS_0           = p.ytrend;
obs_y_FOCUS_25          = p.ytrend;
obs_y_FOCUS_75          = p.ytrend;
obs_y_FOCUS_100         = p.ytrend;
obs_swap_PreDI_3m       = (1 + nrPolicy)^4 - 1;
obs_swap_PreDI_6m       = (1 + nrPolicy)^4 - 1;

obs_y_FOCUS_0_EOQ       = p.ytrend;
obs_y_FOCUS_20_EOQ      = p.ytrend;
obs_y_FOCUS_25_EOQ      = p.ytrend;
obs_y_FOCUS_50_EOQ      = p.ytrend;
obs_y_FOCUS_75_EOQ      = p.ytrend;
obs_y_FOCUS_80_EOQ      = p.ytrend;
obs_y_FOCUS_100_EOQ     = p.ytrend;
obs_unempRate           = 0;
obs_wage                = p.ytrend;

obs_debtOutput_new      = bY;
obs_debtOutput_old      = bY;
obs_tau                 = tau * 100;

%% Planner's objective

%UTIL = 0;
UTIL = -welfare;

%% RETURN

ys = cell2mat({arrayfun(@(x) evalin('caller', x), string(m.endogenous.name))});

% check the validity of the calculations
    %----------------------------------------
if ~utils.error.valid(ys)
    retcode=1;
else
    % push the calculations
    %----------------------
    ss(id)=ys;
end


end


