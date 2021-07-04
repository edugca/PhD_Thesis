%-------------------------------------------------------------
% RBC closed-economy / with and without government
% Reference: Amaral (2020)
% Monetary Policy in a Risky Economy
% Ph.D. Thesis
%-------------------------------------------------------------

endogenous
    c               "$C_t$",
    cR              "$C^R_t$",
    cNR             "$C^{NR}_t$",
    %cStar           "$C^\ast_t$",
    %cRStar          "$C^{R^\ast}_t$",
    %cNRStar         "$C^{{NR}^\ast}_t$",

    @#if conf_appendWelfareEquations
        welfare         "U_t",
        welfareR        "U^R_t",
        welfareNR       "U^{NR}_t",
    @#end

    y               "$Y_t$",
    %yGood           "$Y^{Good}_t$",
    %yBad            "$Y^{Bad}_t$",
    %yGoodExp        "$E_t Y^{Good}_t$",
    %yBadExp         "$E_t Y^{Bad}_t$",
    wp              "$W_t$",
    kk              "$K_t$",
    n               "$N_t$",
    rRN             "$r^{n}_t$",
    rPolicy         "$r_t$",
    rGov            "$r^{Gov}_t$",
    %rGap            "$r^{Gap}_t$",
    %rGapIota        "$r^{Gap,\iota}_t$",
    nrPolicy        "$i_t$",
    nrPolicyExp     "$E_t i_{t+1}$",
    nrGov           "$i^{Gov}_t$",
    iota            "$\overline{\iota}_t$",
    mShock          "$\mathcal{M}_t$",
    mShockExp       "$E_t \mathcal{M}_t$",

    bbeta           "$\tilde{\beta}$",

    @#if endsWith(conf_policyRule, "_priceLevel")
        pLevel      "$P_t$",
    @#end

    Pii             "$\Pi_t$",
    %Pii_YoY         "$\Pi_{t, YoY}$",
    aTilde          "$A_t$",
    aTildeExp       "$E_t A_{t+1}$",
    %recTilde        "$\tilde{\omega}_t$",
    %recTildeExp     "$E_t \tilde{\omega}_t$",
    %def             "$\mathcal{D}^r_t$",
    %defExp          "$E_t \mathcal{D}^r_t$",
    %Psi             "$\Psi_t$",
    %PsiExp          "$E_t \Psi_t$",

    %rGood        "$r^{rn,Good}_t$",
    %rBad         "$r^{rn,Bad}_t$",
    %rGoodOnly    "$r^{rn,GoodOnly}_t$",

    %nrGood        "$i^{Good}_t$",
    %nrBad         "$i^{Bad}_t$",
    %nrGoodOnly    "$i^{GoodOnly}_t$",

    uC              "$U_{C,t}$",
    g               "$G_t$",
    gY              "$\frac{G_t}{Y_t}$"
    z               "$Z_t$",
    b               "$B_t$",
    bLag1           "$B_{t-1}$",
    bY              "$\frac{B_t}{Y_t}$",
    %bd              "$B^d_t$",
    tau             "$\tau_t$",
    tax             "$T_t$",
    polDef          "$\mathcal{D}_t$",
    ddelta          "$\delta_t$",
    muFisLim        "$E_t \mathcal{B}_t$",
    stdFisLim       "$\sigma_t(\mathcal{B}_t)$",
    probDefFisLim   "$Pr_t(\mathcal{B}^{\ast}_{t+1} < B_t)$",
    %drawnFisLim     "$\mathcal{B}^\ast_t$",

    adjCosts    "${AdjCosts}_t$",
    mc          "${MC}_t",

    nEF         "$N^{Eff}_t",
    cREF        "$C^{R,Eff}_t",
    rEF         "r^{Eff}_t",

    nNa         "$N^{N}_t",
    cRNa        "$C^{R,N}_t",
    rNa         "r^{N}_t",

    tauSdw           "$\tau^{sdw}_{t}$",

    %yExp            "$E_t Y_t$",
    %cExp            "$E_t C_t$",
    %cRExp           "$E_t C^R_t$",
    %cNRExp          "$E_t C^{NR}_t$",
    %nExp            "$E_t N_t$",
    %uCExp           "$E_t U_{C,t}$",
    gExp            "$E_t G_{t+1}$",
    %zExp            "$E_t Z_t$",
    %tauExp           "$E_t \tau_{t+1}$",
    %taxExp          "$E_t T_t$",
    ddeltaExp        "$E_t \delta_t$",

    %tauMax          "$\tau^{max}_t$",
    %taxMax          "$T^{max}_t$",
    %wpMax           "$E_t W^{max}_t$",
    %nMax            "$N^{max}_t$",
    %yMax            "$Y^{max}_t$",
    %cMax            "$C^{max}_t$",
    %cRMax           "$C^{R,max}_t$",
    %cNRMax          "$C^{NR,max}_t$",
    %uCMax           "$U_^{max}_{C,t}$",

endogenous
    
    @#if conf_estimateModel || conf_includeObservables
        %obs_pii                 "$\pi^{obs}_t$",
        %obs_nr                  "$i^{obs}_t$",

        @#if conf_modelInterpretation == "Potential_LongRun"
            obs_yPotentialQoQ       "$y^{Potental}_t$",
            obs_lPotentialQoQ       "$n^{Potental}_t$",
            obs_kPotentialQoQ       "$k^{Potental}_t$",
            %obs_y                   "$y^{obs}_t$",
        @#elseif conf_modelInterpretation == "Forecasters"

            @#if ismember('obs_y', conf_observablesList)
                obs_y                   "$y^{obs}_t$",
            @#end
             @#if ismember('obs_c', conf_observablesList)
                obs_c                   "$c^{obs}_t$",
            @#end
            @#if ismember('obs_g', conf_observablesList)
                obs_g                   "$g^{obs}_t$",
            @#end
            @#if ismember('obs_unempRate', conf_observablesList)
                obs_unempRate           "${unemp}^{obs}_t$",
            @#end
            @#if ismember('obs_wage', conf_observablesList)
                obs_wage                   "${\frac{W_t}{P_t}}^{obs}$",
            @#end
            @#if ismember('obs_pii', conf_observablesList)
                obs_pii                  "$\pi^{obs}_t$",
            @#end
            @#if ismember('obs_defLC', conf_observablesList)
                obs_defLC                  "$\mathcal{D}_t$",
            @#end
            @#if ismember('obs_nr', conf_observablesList)
                obs_nr                  "$i_t$",
            @#end
            @#if ismember('obs_swap_PreDI_3m', conf_observablesList)
                obs_swap_PreDI_3m                  "$i_t$",
            @#end
            @#if ismember('obs_swap_PreDI_6m', conf_observablesList)
                obs_swap_PreDI_6m                  "$i_{t+1}$",
            @#end
            @#if ismember('obs_pii_FOCUS_Median', conf_observablesList)
                obs_pii_FOCUS_Median                  "$E_t \pi^{obs}_{t+1}$",
            @#end
            @#if ismember('obs_y_FOCUS_Median_BOQ', conf_observablesList)
                obs_y_FOCUS_Median_BOQ     "$E_tY^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_0', conf_observablesList)
                obs_y_FOCUS_0     "$E_t\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_25', conf_observablesList)
                obs_y_FOCUS_25     "$E_t\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_75', conf_observablesList)
                obs_y_FOCUS_75     "$E_t\tilde{A}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_100', conf_observablesList)
                obs_y_FOCUS_100     "$E_t\tilde{A}^{obs}_t$",
            @#end

            @#if ismember('obs_y_FOCUS_0_EOQ', conf_observablesList)
                obs_y_FOCUS_0_EOQ     "$\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_20_EOQ', conf_observablesList)
                obs_y_FOCUS_20_EOQ     "$\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_25_EOQ', conf_observablesList)
                obs_y_FOCUS_25_EOQ     "$\tilde{\omega}^{obs}_t$",
            @#end
             @#if ismember('obs_y_FOCUS_50_EOQ', conf_observablesList)
                obs_y_FOCUS_50_EOQ     "$\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_75_EOQ', conf_observablesList)
                obs_y_FOCUS_75_EOQ     "$\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_80_EOQ', conf_observablesList)
                obs_y_FOCUS_80_EOQ     "$\tilde{\omega}^{obs}_t$",
            @#end
            @#if ismember('obs_y_FOCUS_100_EOQ', conf_observablesList)
                obs_y_FOCUS_100_EOQ     "$\tilde{A}^{obs}_t$",
            @#end

            @#if ismember('obs_debtOutput_new', conf_observablesList)
                obs_debtOutput_new       "$\left(\frac{B_t}{Y_t}\right)^{obs}$",
            @#end
             @#if ismember('obs_debtOutput_old', conf_observablesList)
                obs_debtOutput_old       "$\left(\frac{B_t}{Y_t}\right)^{obs}$",
            @#end
            @#if ismember('obs_tau', conf_observablesList)
                obs_tau                  "$\tau^{obs}$",
            @#end

        @#else
            obs_y                   "$y^{obs}_t$",
            %obs_c                  "$c^{obs}_t$",
        @#end

        %obs_pii_FOCUS_median,
        %obs_y_FOCUS_median,
    @#end

exogenous
    epsA            "$\varepsilon^A_t$",
    epsM            "$\varepsilon^M_t$",

    @#if ~ismember('sigmaBeta', conf_shocksTurnedOff)
        epsBeta         "$\varepsilon^\beta_t$",
    @#end

    @#if conf_shockOnEffectiveProductivity ~= "None"
        epsPsi          "$\varepsilon^{\Psi}_t$",
    @#end

    @#if conf_prodProcessType == "Dual_LogNormal" || conf_prodProcessType == "MarkovSwitching"
        epsR           "$\varepsilon^R_t$",
    @#end

    @#if conf_defRProcessType == "LogNormal"
        epsD           "$\varepsilon^D_t$",
    @#end

    @#if conf_kVariable == "Exogenous"
        epsK        "$\varepsilon^K_t$",
    @#end

    %@#if conf_hasGovernment
        epsG        "$\varepsilon^G_t$",

        @#if ~ismember('sigmaTau', conf_shocksTurnedOff)
            epsTau      "$\varepsilon^\tau_t$",
        @#end

        %@#if conf_govBondsAreRiskFree == false
            %epsPolDef       "$\varepsilon^\mathcal{D}_t$",
            %epsFL           "$\varepsilon^\mathcal{B}_t$",
        %@#end
    %@#end

parameters

    gammaTau        "$\gamma^{\tau}$",
    phi             "$\phi^{\pi}$",

    elbRate         "\overline{i^{ELB}}",

    debtWedge,

    beta            "$\beta$",
    sigma           "$\sigma$",
    eta             "$\eta$",
    chi             "$\chi$",
    alphaG          "$\alpha_G$"
    gammaGPSI       "$\gamma_{G\Psi}$"

    phiNr           "$\phi^i$",
    phi_Y           "$\phi^Y$",
    phi_dY          "$\phi^{dY}$",
    phi_Exp         "$\phi^{Exp}$",
    pLevelBar       "$\overline{P}$",
    piiBar          "$\overline{\pi}$",
    
    rhoA            "$\rho^A$",
    rhoR            "$\rho^R$",
    rhoD            "$\rho^D$",
    rhoM            "$\rho^{\mathcal{M}}$",
    rhoBeta         "$\rho^\beta$",
    rhoPolDef       "$\rho^{\mathcal{D}}$",

    rhoAG            "$\rho^{AG}$",

    sigmaA          "$\sigma^A$",
    sigmaR          "$\sigma^R$",
    sigmaD          "$\sigma^D$",
    sigmaM          "$\sigma^M$",
    sigmaK          "$\sigma^K$",
    sigmaPsi        "$\sigma^\Psi$",
    sigmaBeta       "$\sigma^\beta$",
    sigmaPolDef     "$\sigma^{\mathcal{D}}$",
    
    %aTildeSS        "$\overline{A}$",
    aTildeBar       "$\overline{A}$",
    recTildeSS      "$\overline{\omega}$",
    defSS           "$\overline{\mathcal{D}}$",
    mShockSS        "$\overline{\mathcal{M}}$",

    ySS             "$\overline{Y}$",

    kBar            "$\overline{K}$",

    fracNR          "$\omega^{NR}$",

    ytrend          "$\overline{dY/dt}$",

    unempSS         "$\overline{U}$",

    tfpLoss,        "$A^{loss}$",

    %%%%%%%%%%%%%%%% GOVERNMENT
    rhoTau          "$\rho^{\tau}$",
    sigmaTau        "$\sigma^\tau$",
    rhoGG           "$\rho^{GG}$",
    rhoGY           "$\rho^{GY}$",
    sigmaG          "$\sigma^{G}$",

    bbetaSS         "$\overline{\beta}$",
    yBar            "$\overline{Y}$", % steady-state y at REGIME 1
    gYSS            "$\overline{G/Y}$",
    bYSS            "$\overline{B/Y}$",
    zYSS            "$\overline{Z/Y}$",
    %tauSS           "$\overline{\tau}$",
    %tauMaxSS        "$\overline{\tau^{max}}$",
    polDefSS        "$\overline{\mathcal{D}}$",

    tauBar          "$\overline{\tau}$",
    tauMaxBar       "$\overline{\tau^{max}}$",
    gYBar           "$\overline{G/Y}$",
    bYBar           "$\overline{B/Y}$",
    zYBar           "$\overline{Z/Y}$",

    deltaBar        "$\overline{\delta}$",

    polDefWedge     "$\overline{\mathcal{D}^{\ast}} - \overline{\mathcal{D}}$",


    thetaElast      "$\theta^{elast}$",
    phiC            "$\phi^C$",
    labSub          "${labSub}$",

    %%%%%%%%%%%%%%%%

% Fiscal limits parameters
%---------------------------------------------------------------------------
parameters
    muFisLim_Intercept,
    muFisLim_aTilde,
    muFisLim_g,
    %muFisLim_bbeta,

    stdFisLim_Intercept,
    stdFisLim_aTilde,
    stdFisLim_g,
    %stdFisLim_bbeta,

    probDefFisLim_Param_0,
    probDefFisLim_Param_B,
    probDefFisLim_Param_A,
    probDefFisLim_Param_G,
    %probDefFisLim_Param_Beta


% measurement errors are declared as parameters: The corresponding variables
% have to be declared as observables!!!
%---------------------------------------------------------------------------
parameters
    @#if conf_estimateModel

        @#if conf_measurementErrors
            @#if conf_modelInterpretation == "Potential_LongRun"
                stderr_obs_pii,
                stderr_obs_yPotentialQoQ,
                stderr_obs_lPotentialQoQ,
                @#if conf_kVariable == "Exogeonus"
                    stderr_obs_kPotentialQoQ,
                @#end
            @#else conf_modelInterpretation == "Forecasters"
                @#if ismember('obs_y', conf_measurementErrorsList)
                    stderr_obs_y            "$\sigma^{me,Y}$",
                @#end
                @#if ismember('obs_c', conf_measurementErrorsList)
                    stderr_obs_c,
                @#end
                @#if ismember('obs_g', conf_measurementErrorsList)
                    stderr_obs_g,
                @#end
                @#if ismember('obs_unempRate', conf_measurementErrorsList)
                    stderr_obs_unempRate        "$\sigma^{me,N}$",
                @#end
                @#if ismember('obs_wage', conf_measurementErrorsList)
                    stderr_obs_wage             "$\sigma^{me,W}$",
                @#end
                @#if ismember('obs_pii', conf_measurementErrorsList)
                    stderr_obs_pii,
                @#end
                @#if ismember('obs_defLC', conf_measurementErrorsList)
                    stderr_obs_defLC            "$\sigma^{me,\mathcal{D}}$",
                @#end
                @#if ismember('obs_nr', conf_measurementErrorsList)
                    stderr_obs_nr,
                @#end
                @#if ismember('obs_swap_PreDI_3m', conf_measurementErrorsList)
                    stderr_obs_swap_PreDI_3m,
                @#end
                @#if ismember('obs_swap_PreDI_6m', conf_measurementErrorsList)
                    stderr_obs_swap_PreDI_6m,
                @#end
                @#if ismember('obs_pii_FOCUS_Median', conf_measurementErrorsList)
                    stderr_obs_pii_FOCUS_Median,
                @#end
                @#if ismember('stderr_obs_y_FOCUS_Median_BOQ', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_Median_BOQ,
                @#end
                @#if ismember('obs_y_FOCUS_0', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_0,
                @#end
                @#if ismember('obs_y_FOCUS_25', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_25,
                @#end
                @#if ismember('obs_y_FOCUS_75', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_75,
                @#end
                @#if ismember('obs_y_FOCUS_100', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_100,
                @#end

                @#if ismember('obs_y_FOCUS_20_EOQ', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_20_EOQ "$\sigma^\tilde{\omega}_{\text{error}}$",
                @#end
                @#if ismember('obs_y_FOCUS_80_EOQ', conf_measurementErrorsList)
                    stderr_obs_y_FOCUS_80_EOQ "$\sigma^\tilde{A}_{\text{error}}$",
                @#end

            @#end
        @#end
    @#end

@#if conf_occBinConstraint == "FiscalLimit"
    %parameters r1fisLim_tp_1_2 r1fisLim_tp_2_1
    %parameters r2taxLim_tp_1_2 r2taxLim_tp_2_1

    parameters(r1fisLim,2)
        bindFis "$\mathbbm{1}_{\text{fis limit}}$"

    parameters(r2taxLim,2)
        bindTax "$\mathbbm{1}_{\text{tax limit}}$"

@#else
    parameters
        bindTax
        bindFis
@#end

@#if conf_effectiveLowerBound
     parameters(r3elbLim,2)
        bindELB "$\mathbbm{1}_{\text{elb limit}}$",
@#else
    parameters
        bindELB
@#end

model

    % Auxiliary
    # aTildeSS = aTildeBar;
    # aTildeExpSS = aTildeBar;
    # bSS   = bYBar * yBar;
    # gSS   = gYBar * yBar;
    # zSS   = zYBar * yBar;
    # tauSS = tauBar;
    # tauSdwSS = tauBar;
    # tauExpSS = tauBar;
    # tauMaxSS = tauMaxBar;

    % Regime switching and cccasionally binding constraints

        % Endogenous switching probabilities
        @#if conf_occBinConstraint == "FiscalLimit"
            % Default probabilities
            ! r1fisLim_tp_1_2 = 1/(1 + exp( probDefFisLim_Param_0 + probDefFisLim_Param_B*(b - muFisLim) + probDefFisLim_Param_A*(aTildeExp-aTildeExpSS) + probDefFisLim_Param_G*(gExp-gSS) )) ;
            %! r1fisLim_tp_2_1 = 1 - 1/(1 + exp( probDefFisLim_Param_0 + probDefFisLim_Param_B*(b - muFisLim) + probDefFisLim_Param_A*(aTildeExp-aTildeExpSS) + probDefFisLim_Param_G*(gExp-gSS) )) ;
            ! r1fisLim_tp_2_1 = 1 ;

            ! r2taxLim_tp_1_2 = (tauSdw > 0.5) ; %0.5 is the maximum tax rate
            ! r2taxLim_tp_2_1 = (tauSdw < 0.5) ;

            % Check whether tau is larger than tauMax
            %? tau < tauMaxSS ;

            % Check whether the government has defaulted
            %? bLag1 < drawnFisLim ;
        @#end

        @#if conf_effectiveLowerBound
            ! r3elbLim_tp_1_2 = (nrPolicyExp < elbRate) ;
            ! r3elbLim_tp_2_1 = (nrPolicyExp > elbRate) ;
        @#end

        @#if conf_hasGovernment
            bLag1   = b(-1);

            @#if conf_debtLevel == "Zero"
                @#if ismember('sigmaTau', conf_shocksTurnedOff)
                    %log(tau) = (1-bindTax)*(log(tauSS) + rhoTau*log(tau(-1)/tauSS) ) + bindTax*(log(tauMaxSS));
                    log(tau) = (1-bindTax)*( log(tauSdw(-1)) ) + bindTax*(log(tauMaxSS));
                @#else
                    %log(tau) = (1-bindTax)*(log(tauSS) + rhoTau*log(tau(-1)/tauSS) + sigmaTau*epsTau) + bindTax*(log(tauMaxSS));
                    log(tau) = (1-bindTax)*( log(tauSdw(-1)) ) + bindTax*(log(tauMaxSS));
                @#end
                log(tauSdw) = log(tauSdwSS) + rhoTau*log(tauSdw(-1)/tauSdwSS) ;
            @#else
                @#if ismember('sigmaTau', conf_shocksTurnedOff)
                    @#if conf_taxRule == "Relative"
                        %log(tau) = (1-bindTax)*(log(tauSS) + rhoTau*log(tau(-1)/tauSS) + gammaTau*log(b(-1)/y(-1) * steady_state(y)/steady_state(b)) ) + bindTax*(log(tauMaxSS));
                        log(tau) = (1-bindTax)*( log(tauSdw(-1)) ) + bindTax*(log(tauMaxSS));
                    @#elseif conf_taxRule == "Absolute"
                        %log(tau) = (1-bindTax)*(log(tauSS) + rhoTau*log(tau(-1)/tauSS) + gammaTau*log(b(-1)/steady_state(b)) ) + bindTax*(log(tauMaxSS));
                        log(tau) = (1-bindTax)*( log(tauSdw(-1)) ) + bindTax*(log(tauMaxSS));
                    @#end
                @#else
                    @#if conf_taxRule == "Relative"
                        %log(tau) = (1-bindTax)*(log(tauSS) + rhoTau*log(tau(-1)/tauSS) + gammaTau*log(b(-1)/y(-1) * steady_state(y)/steady_state(b)) + sigmaTau*epsTau) + bindTax*(log(tauMaxSS));
                        log(tau) = (1-bindTax)*( log(tauSdw(-1)) + sigmaTau*epsTau ) + bindTax*(log(tauMaxSS));
                    @#elseif conf_taxRule == "Absolute"
                        %log(tau) = (1-bindTax)*(log(tauSS) + rhoTau*log(tau(-1)/tauSS) + gammaTau*log(b(-1)/steady_state(b)) + sigmaTau*epsTau) + bindTax*(log(tauMaxSS));
                        log(tau) = (1-bindTax)*( log(tauSdw(-1)) + sigmaTau*epsTau ) + bindTax*(log(tauMaxSS));
                    @#end
                @#end
            
                @#if conf_taxRule == "Relative"
                    %log(tauSdw) = log(tauSdwSS) + rhoTau*log(tauSdw(-1)/tauSdwSS) + gammaTau*log(b/y * steady_state(y)/steady_state(b)) ;
                    log(tauSdw) = log(tauSdwSS) + rhoTau*log(tauSdw(-1)/tauSdwSS) + gammaTau*log(b/y * ySS/bSS) ;
                @#elseif conf_taxRule == "Absolute"
                    log(tauSdw) = log(tauSdwSS) + rhoTau*log(tauSdw(-1)/tauSdwSS) + gammaTau*log(b/bSS) ;
                @#end
            @#end

             @#if conf_govBondsAreRiskFree
                % Set fiscal limits to infinity
                muFisLim = 0 ;
                stdFisLim = 0 ;
                probDefFisLim = 0;
                %drawnFisLim = 0;
                ddelta = 0 ;
                ddeltaExp = 0;
            @#else
                % Default?
                %muFisLim = fFiscalLimit_mu(aTilde,recTilde,def,g,bbeta);
                %stdFisLim = fFiscalLimit_std(aTilde,recTilde,def,g,bbeta);
                %@#include "rise_fiscalLimit.rs"

                %[name='Expected fiscal limit mean']
                muFisLim  = muFisLim_Intercept + muFisLim_aTilde*(aTildeExp-aTildeExpSS) + muFisLim_g*(gExp-gSS) ;
                %[name='Expected fiscal limit standard deviation']
                stdFisLim = stdFisLim_Intercept + stdFisLim_aTilde*(aTildeExp-aTildeExpSS) + stdFisLim_g*(gExp-gSS) ;
                %[name='Expected default probability']
                probDefFisLim = 1/(1 + exp( probDefFisLim_Param_0 + probDefFisLim_Param_B*(b - muFisLim) + probDefFisLim_Param_A*(aTildeExp-aTildeExpSS) + probDefFisLim_Param_G*(gExp-gSS) )) ;

                %drawnFisLim = muFisLim + stdFisLim*epsFL;
                %probDefFisLim = fFiscalLimit_defProb(bLag1, aTilde-steady_state(aTilde), recTilde - steady_state(recTilde), def - steady_state(def), g - steady_state(g), bbeta - steady_state(bbeta));
                %probDefFisLim = normcdf(bLag1, muFisLim, stdFisLim) ;
                %probDefFisLim = sigmaPolDef*epsPolDef;
                
                %@#if strcmp(conf_approximateDefaultProb, 'logistic')
                %    probDefFisLim = 1/(1 + exp( probDefFisLim_Param_0 + probDefFisLim_Param_B*(bLag1 - muFisLim) + probDefFisLim_Param_A*(aTilde-aTildeSS) + probDefFisLim_Param_G*(g-gSS) + sigmaPolDef*epsPolDef ));
                %@#elseif strcmp(conf_approximateDefaultProb, 'linear')
                %    probDefFisLim = probDefFisLim_Param_0 + probDefFisLim_Param_B*(bLag1 - muFisLim) + probDefFisLim_Param_A*(aTilde-aTildeSS) + probDefFisLim_Param_G*(g-gSS) + sigmaPolDef*epsPolDef;
                %@#end
                
                ddelta       = bindFis * deltaBar ;
                %deltaExp    = polDef * deltaBar ;
                ddeltaExp    = ddelta(+1);

            @#end

        @#else
            bLag1           = 0;
            tau             = 0;
            tauSdw          = 0;
            %tauExp          = 0;
            ddelta           = 0;
            ddeltaExp        = 0;
            muFisLim        = 0;
            stdFisLim       = 0;
            %probDefFisLim   = 0;
            %drawnFisLim     = 0;
        @#end

	% Exogenous processes

        %[name='aTilde process']
        %TFP loss just in the first period that the fiscal constraint binds: (1-polDef(-1))
        log(aTilde/aTildeSS)    = rhoA*log(aTilde(-1)/aTildeSS) + (1 - bindFis)*(sigmaA*epsA + rhoAG*epsG) - bindFis*(1-polDef(-1))*tfpLoss;
        log(aTildeExp/aTildeExpSS) = rhoA*log(aTilde/aTildeSS) ;

        %[name='Monetary shock process']
        log(mShock/mShockSS)      = rhoM*log(mShock(-1)/mShockSS) + sigmaM*epsM ;
        %[name='Expected monetary shock process']
        log(mShockExp/mShockSS)   = rhoM*log(mShock/mShockSS) ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GOVERNMENT EQUATIONS
        
        @#if conf_hasGovernment
            @#if conf_govExpenses == "Zero"
                g       = rhoGG*g(-1) + rhoGY*log(y(-1)/ySS) + sigmaG*epsG;
                gExp    =  rhoGG*g + rhoGY*log(y/ySS);
            @#else
                %[name='g process']
                log(g/gSS)      = rhoGG*log(g(-1)/gSS) + rhoGY*log(y(-1)/ySS) + sigmaG*epsG;
                log(gExp/gSS)   = rhoGG*log(g/gSS) + rhoGY*log(y/ySS);
            @#end
            gY              = g/y;

            %[name='Tranfers']
            z       = zSS;

            @#if conf_govAccumDebt
                %[name='Debt']
                %b       = (1 + nrGov) * ( (1 - bindFis * deltaBar) * b(-1)/Pii + g + z + debtWedge - tax);
                % bindFis*(1-polDef(-1))*deltaBar
                
                b       =   (1-bindFis) * ( (1 + nrGov) * ( (1 - 0) * b(-1)/Pii + g + z + debtWedge - tax) )
                            + bindFis * ( (1 + nrGov) * ( (1 - deltaBar) * b(-1)/Pii + g + z + debtWedge - tax) );
            @#else
                %[name='Debt']
                b       = 0;
            @#end

            %bd      = b(-1);
            bY      = b/y;
        @#else
            g       = 0;
            gY      = 0;
            %gExp    = 0;
            
            z       = 0;
            %zExp    = 0;
        
            b       = 0;
            %bd      = 0;
            bY      = 0;
        @#end
        
        
    % CPOs of Households
        
        @#if conf_appendWelfareEquations
            %[name='Welfare of the Ricardian household']
            welfareR    = 1/(1-sigma)*(cR + alphaG*g - eta*n^(1+chi)/(1+chi))^(1-sigma) + bbeta*welfareR(+1);
            %[name='Welfare of the non-Ricardian household']
            welfareNR   = 1/(1-sigma)*(cNR + alphaG*g - eta*n^(1+chi)/(1+chi))^(1-sigma) + bbeta*welfareNR(+1);
            %[name='Welfare']
            welfare     = (1-fracNR)*welfareR + fracNR*welfareNR;
        @#end

        %[name='Marginal utility of consumption of the Ricardian household']
        uC = ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) )^(-sigma);
        
    %%%%%%%%%%%%%%% FLEXIBLE/STICKY PRICES
        %[name='Taxes']
        tax = tau*(n*wp) + tau*(y - n*wp - adjCosts);

        %[name='Adjustment costs']
        adjCosts = phiC/2 * (Pii/(steady_state(Pii)) - 1)^2 * y;

        %[name='Marginal cost']
        mc = wp / (aTilde*kk);

        %[name='Non-linear NK Phillips curve']
        (1-thetaElast) + thetaElast*mc - phiC*(Pii/(steady_state(Pii)) - 1)*Pii/(steady_state(Pii)) + phiC * bbeta*uC(+1)/uC * (Pii(+1)/(steady_state(Pii)) - 1)*Pii(+1)/(steady_state(Pii))*y(+1)/y = 0;

        % Equilibrium condition
        %[name='Resources constraint']
        c       = y - g - adjCosts;

        %[name='Production function']
        y           = aTilde * kk * n;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %tax = tau * kk * ((1-tau) * kk/eta * aTilde)^(1/chi) * aTilde;
        %% CPOs of Firms
        %%[name='Labor demand']
        %n   = ( (1-tau) * kk/eta * aTilde )^(1/chi);
        %% Equilibrium condition
        %%[name='Resources constraint']
        %c       = y - g;
        %%cStar   = y - (1-alphaG)*g;
        %%[name='Production function']
        %y           = kk * ( (1-tau) * kk/eta * aTilde )^(1/chi) * aTilde;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %[name='Fisher equation']
        %(1+nr) = (1+rRN)*Pii(+1);
        @#if startsWith(conf_policyRule, "rRF_")
            (1+nrPolicy) = (1+rRN)*Pii(+1);
        @#else
            (1+nrPolicy) = (1+rPolicy)*Pii(+1);
        @#end
        %(1+nrPolicy) = (1-polDef(+1))*(1+rPolicy)*Pii(+1) + polDef(+1)*(1+rPolicy)*(1-ddelta(+1))*Pii(+1);
        
        %[name='Nominal interest rate of government bond']
        (1+nrGov) = bindELB*(1 + elbRate) + (1-bindELB)*(1+rGov)*Pii(+1);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CPOs of Households
    
        %[name='Labor supply']
        wp = (1/(1-tau)) * eta * n^chi ;

        %[name='Household Euler (risk-free)']
        1 / (1 + rRN) = bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;

        @#if conf_govBondsAreRiskFree
            %[name='Household Euler (policy is risk-free)']
            1 / (1 + rPolicy) =  bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;

            %[name='Household Euler (gov bond)']
            rGov = rPolicy;
        @#elseif startsWith(conf_policyRule, "rRF_") || startsWith(conf_policyRule, "rRNwithRiskFree_")
            %[name='Household Euler (policy is risk-free)']
            1 / (1 + rPolicy) =  bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;

            %[name='Household Euler (gov bond)']
            1 / (1 + rGov) = (1 - probDefFisLim) * bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 
                             +  ( probDefFisLim*(1-deltaBar) ) * bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;
        @#else
            %[name='Household Euler (policy is risky)']
            1 / (1 + rPolicy) = (1 - probDefFisLim) * bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 
                             +  ( probDefFisLim*(1-deltaBar) ) * bbeta * ( ( cR(+1) + alphaG*g(+1) - eta*n(+1)^(1+chi)/(1+chi) ) / ( cR + alphaG*g - eta*n^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;
            
            %[name='Household Euler (gov bond)']
            rGov = rPolicy;
        @#end

    % Natural real interest rate
        %[name='Employment (natural)']
        nNa = ( (thetaElast - 1)/thetaElast * (1-tau)/eta * kk * aTilde  )^(1/chi) ;

        %[name='Ricardian consumption (natural)']
        cRNa = 1/(1-fracNR) * ( ( (thetaElast - 1)/thetaElast * (1-tau)/eta )^(1/chi) * (kk*aTilde)^(1+1/chi) - g - fracNR*(( (thetaElast - 1)/thetaElast * 1/eta)^(1/chi) * ((1-tau)*kk*aTilde)^(1+1/chi) + z ) ) ;

        %[name='Household Euler (risk-free and natural)']
        1 / (1 + rNa) = bbeta * ( ( cRNa(+1) + alphaG*g(+1) - eta*nNa(+1)^(1+chi)/(1+chi) ) / ( cRNa + alphaG*g - eta*nNa^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;

    % Efficient real interest rate
        %[name='Employment (efficient)']
        nEF = ( (1-tau)/eta * kk * aTilde  )^(1/chi) ;

        %[name='Ricardian consumption (efficient)']
        cREF = 1/(1-fracNR) * ( ( (1-tau)/eta )^(1/chi) * (kk*aTilde)^(1+1/chi) - g - fracNR*((1/eta)^(1/chi) * ((1-tau)*kk*aTilde)^(1+1/chi) + z ) ) ;

        %[name='Household Euler (risk-free and efficient)']
        1 / (1 + rEF) = bbeta * ( ( cREF(+1) + alphaG*g(+1) - eta*nEF(+1)^(1+chi)/(1+chi) ) / ( cREF + alphaG*g - eta*nEF^(1+chi)/(1+chi) ) )^(-sigma) * 1 ;

    % Policy default probability
        @#if conf_hasGovernment == true
            @#if conf_govBondsAreRiskFree == false
                polDef = bindFis;
            @#else
                polDef = 0;
            @#end
        @#else
            polDef = 0;
        @#end

        % Future tax rates
        %yExp     = kk * ( (1-tauExp) * kk/eta * PsiExp )^(1/chi) * PsiExp;
        %cExp     = yExp - gExp;
        %cNRExp   = (1-tauExp) * wp * nExp + zExp + debtWedge;
        %cExp     = fracNR * cNRExp + (1-fracNR) * cRExp ; 

        %uCExp    = ( cExp + alphaG*gExp - eta*nExp^(1+chi)/(1+chi) )^(-sigma);
        %uCExp    = ( cRExp + alphaG*gExp - eta*nExp^(1+chi)/(1+chi) )^(-sigma);

        %nExp     = ( (1-tauExp) * kk/eta * PsiExp )^(1/chi);
        %taxExp   = tau * (1-tau)^(1/chi) * (1/eta)^(1/chi) * (kk * PsiExp)^(1+1/chi);

    % Capital
        %[name='Auxiliary variable: Capital accumulation']
        @#if conf_kVariable == "Fixed"
            @#if conf_hasGovernment
                @#if conf_govExpenses == "Zero"
                    kk = kBar * (1 + g - gSS)^gammaGPSI;
                @#else
                    kk = kBar * (g/gSS)^gammaGPSI;
                @#end
            @#else
                kk = kBar ;
            @#end
        @#else
            @#if conf_hasGovernment
                @#if conf_govExpenses == "Zero"
                    kk = kk(-1) * (1 + g - gSS)^gammaGPSI + sigmaK * epsK;
                @#else
                    kk = kk(-1) * (g/gSS)^gammaGPSI + sigmaK * epsK;
                @#end
            @#else
                kk = kk(-1) + sigmaK * epsK;
            @#end
        @#end

    % Non-Ricardian households
        %[name='Consumption of non-Ricardian households']
        cNR     = (1-tau) * wp * n + z + debtWedge;
        %[name='Aggregate consumption']
        c       = fracNR * cNR + (1-fracNR) * cR ;

    % Policy rule

        @#if ~ismember('sigmaBeta', conf_shocksTurnedOff)
            %[name='Preferences shock']
            log(bbeta/beta) = rhoBeta*log(bbeta(-1)/beta) + sigmaBeta*epsBeta;
        @#else
            %[name='Preferences shock']
            bbeta = beta;
        @#end
    
        @#if endsWith(conf_policyRule, "_priceLevel")
            %[name='Price level']
            pLevel = pLevel(-1) * Pii;
        @#end

        %[name='YoY Inflation']
        %Pii_YoY = Pii/Pii(-4) ;

        %[name='rGap']
        %rGap        = rPolicy - rRN ;
        %rGapIota    = iota - rRN ;

        @#if endsWith(conf_policyRule, "_inflation")
            %[name='Central bank rule']
            (1 + nrPolicy) =
                (1-bindELB) * (1 + nrPolicy(-1))^phiNr * ( (1 + iota) * (steady_state(Pii)) * ( Pii/(steady_state(Pii)) )^phi * ( y(-1)/steady_state(y) )^phi_Y * ( y(-1)/y(-2) )^phi_dY )^(1-phiNr) * mShock
                + bindELB * (1 + elbRate)
                ;
        
            %[name='ELB: expected policy rule']
            (1 + nrPolicyExp) = (1 + nrPolicy)^phiNr * ( (1 + iota(+1)) * (steady_state(Pii)) * ( Pii(+1)/(steady_state(Pii)) )^phi * ( y/steady_state(y) )^phi_Y * ( y/y(-1) )^phi_dY )^(1-phiNr) * mShockExp;
        @#elseif endsWith(conf_policyRule, "_priceLevel")
            %[name='Central bank rule']
            (1 + nrPolicy) = 
                (1-bindELB) * (1 + nrPolicy(-1))^phiNr * ( (1 + iota) * (steady_state(Pii)) * ( pLevel/pLevelBar )^phi * ( y(-1)/steady_state(y) )^phi_Y * ( y(-1)/y(-2) )^phi_dY )^(1-phiNr) * mShock
                + bindELB *(1 + elbRate)
                ;

             %[name='ELB: expected policy rule']
             (1 + nrPolicyExp) = (1 + nrPolicy)^phiNr * ( (1 + iota(+1)) * (steady_state(Pii)) * ( pLevel(+1)/pLevelBar )^phi * ( y/steady_state(y) )^phi_Y * ( y/y(-1) )^phi_dY )^(1-phiNr) * mShockExp;
        @#end
        
        %[name='Central bank rule: intercept']
        @#if conf_govBondsAreRiskFree
            @#if startsWith(conf_policyRule, "fixedIntercept_")
                iota = -1 + 1/beta ;
            @#elseif startsWith(conf_policyRule, "rRN_")
                iota = rRN ;
            @#elseif startsWith(conf_policyRule, "rRF_")
                iota = rRN ;
            @#elseif startsWith(conf_policyRule, "rEF_")
                iota = rEF ;
            @#elseif startsWith(conf_policyRule, "rNa_")
                iota = rNa ;
            @#elseif startsWith(conf_policyRule, "polDefAdjusted_")
                iota = rPolicy ;
            @#elseif startsWith(conf_policyRule, "rRNwithRiskFree_")
                iota = rGov ;
            @#end
        @#else
            @#if startsWith(conf_policyRule, "fixedIntercept_")
                iota = steady_state(rPolicy) ;
            @#elseif startsWith(conf_policyRule, "rRN_")
                iota = rRN + (steady_state(rGov) - steady_state(rRN)) ;
            @#elseif startsWith(conf_policyRule, "rRF_")
                iota = rRN ;
            @#elseif startsWith(conf_policyRule, "rEF_")
                iota = rEF + (steady_state(rGov) - steady_state(rRN)) ;
            @#elseif startsWith(conf_policyRule, "rNa_")
                iota = rNa + (steady_state(rGov) - steady_state(rRN)) ;
            @#elseif startsWith(conf_policyRule, "polDefAdjusted_")
                iota = rPolicy ;
             @#elseif startsWith(conf_policyRule, "rRNwithRiskFree_")
                iota = rGov - (steady_state(rGov) - steady_state(rRN)) ;
            @#end
        @#end

    % measurement equations
    @#if conf_estimateModel || conf_includeObservables
        
        @#if conf_modelInterpretation == "Potential_LongRun"
            % Potential model

            obs_yPotentialQoQ           = log(aTilde/aTilde(-1));
            obs_lPotentialQoQ           = log(n/n(-1));
            @#if conf_kVariable == "Exogenous"
                obs_kPotentialQoQ       = log(kk/kk(-1));
            @#end

            %obs_pii_FOCUS_median    = Pii(+1) - 1;
            %obs_y_FOCUS_median      = log(y(+1)/y(-3));

        @#elseif conf_modelInterpretation == "Forecasters"
            @#if ismember('obs_y', conf_observablesList)
                obs_y                       = log(y/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_c', conf_observablesList)
                obs_c                       = log(c/c(-1)) + ytrend;
            @#end
            @#if ismember('obs_g', conf_observablesList)
                @#if conf_govExpenses == "Zero"
                    obs_g                       = 0;
                @#else
                    obs_g                       = log(g/g(-1)) + ytrend;
                @#end
            @#end
            @#if ismember('obs_unempRate', conf_observablesList)
                obs_unempRate  = n/steady_state(n) - 1 ;
            @#end
            @#if ismember('obs_wage', conf_observablesList)
                obs_wage                       = log(wp/wp(-1)) + ytrend;
            @#end
            @#if ismember('obs_pii', conf_observablesList)
                obs_pii                     = log(Pii);
            @#end
            @#if ismember('obs_defLC', conf_observablesList)
                obs_defLC                     = (polDef-polDef(-1)) + polDefWedge;
            @#end
            @#if ismember('obs_nr', conf_observablesList)
                obs_nr                     = (1 + nrPolicy)^4 - 1;
            @#end
            @#if ismember('obs_swap_PreDI_3m', conf_observablesList)
                obs_swap_PreDI_3m          = (1 + nrPolicy)^4 - 1;
            @#end
            @#if ismember('obs_swap_PreDI_6m', conf_observablesList)
                obs_swap_PreDI_6m          = ( (1+nrPolicy)*(1 + nrPolicy(+1)) )^2 - 1;
            @#end
             @#if ismember('obs_pii_FOCUS_Median', conf_observablesList)
                obs_pii_FOCUS_Median/100   =  log(Pii(+1));
            @#end
            @#if ismember('obs_y_FOCUS_Median_BOQ', conf_observablesList)
                obs_y_FOCUS_Median_BOQ   =  log(yExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_0', conf_observablesList)
                obs_y_FOCUS_0   = log(yBadExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_25', conf_observablesList)
                obs_y_FOCUS_25   = log(yBadExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_75', conf_observablesList)
                obs_y_FOCUS_75   = log(yGoodExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_100', conf_observablesList)
                obs_y_FOCUS_100   = log(yGoodExp/y(-1)) + ytrend;
            @#end

            @#if ismember('obs_y_FOCUS_0_EOQ', conf_observablesList)
                obs_y_FOCUS_0_EOQ   = log(yBad/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_20_EOQ', conf_observablesList)
                obs_y_FOCUS_20_EOQ   = log(yBadExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_25_EOQ', conf_observablesList)
                obs_y_FOCUS_25_EOQ   = log(yBadExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_50_EOQ', conf_observablesList)
                obs_y_FOCUS_50_EOQ   = log(y/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_75_EOQ', conf_observablesList)
                obs_y_FOCUS_75_EOQ   = log(yGoodExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_80_EOQ', conf_observablesList)
                obs_y_FOCUS_80_EOQ   = log(yGoodExp/y(-1)) + ytrend;
            @#end
            @#if ismember('obs_y_FOCUS_100_EOQ', conf_observablesList)
                obs_y_FOCUS_100_EOQ   = log(yGood/y(-1)) + ytrend;
            @#end

            @#if ismember('obs_debtOutput_new', conf_observablesList)
                obs_debtOutput_new/100 = bY;
            @#end
             @#if ismember('obs_debtOutput_old', conf_observablesList)
                obs_debtOutput_old/100 = bY;
            @#end
            @#if ismember('obs_tau', conf_observablesList)
                obs_tau/100 = tau;
            @#end

        @#end

        %obs_yPotentialQoQ   = log(aTilde/aTilde(-1));
        %obs_y               = log(y/y(-1));
		%obs_c               = log(c/c(-1));

        %obs_pii             = Pii - 1;
        %obs_nr              = (1 + nrPolicy)^4 - 1;
    @#end

observables

    @#if conf_estimateModel || conf_includeObservables
        @#if conf_modelInterpretation == "Potential_LongRun"
            obs_yPotentialQoQ
            obs_lPotentialQoQ
            @#if conf_kVariable == "Exogenous"
                obs_kPotentialQoQ
            @#end
        @#elseif conf_modelInterpretation == "Forecasters"
            @#if ismember('obs_y', conf_observablesList)
                obs_y,
            @#end
            @#if ismember('obs_c', conf_observablesList)
                obs_c,
            @#end
            @#if ismember('obs_g', conf_observablesList)
                obs_g,
            @#end
            @#if ismember('obs_unempRate', conf_observablesList)
                obs_unempRate,
            @#end
            @#if ismember('obs_wage', conf_observablesList)
                obs_wage,
            @#end
            @#if ismember('obs_pii', conf_observablesList)
                obs_pii,
            @#end
            @#if ismember('obs_defLC', conf_observablesList)
                obs_defLC,
            @#end
            @#if ismember('obs_nr', conf_observablesList)
                obs_nr,
            @#end
            @#if ismember('obs_swap_PreDI_3m', conf_observablesList)
                obs_swap_PreDI_3m,
            @#end
            @#if ismember('obs_swap_PreDI_6m', conf_observablesList)
                obs_swap_PreDI_6m,
            @#end
            @#if ismember('obs_pii_FOCUS_Median', conf_observablesList)
                obs_pii_FOCUS_Median,
            @#end
            @#if ismember('obs_y_FOCUS_Median_BOQ', conf_observablesList)
                obs_y_FOCUS_Median_BOQ,
            @#end
            @#if ismember('obs_y_FOCUS_0', conf_observablesList)
                obs_y_FOCUS_0,
            @#end
            @#if ismember('obs_y_FOCUS_25', conf_observablesList)
                obs_y_FOCUS_25,
            @#end
            @#if ismember('obs_y_FOCUS_75', conf_observablesList)
                obs_y_FOCUS_75,
            @#end
            @#if ismember('obs_y_FOCUS_100', conf_observablesList)
                obs_y_FOCUS_100,
            @#end
            
            %%%%%%%%%%%%%%% EOQ
            @#if ismember('obs_y_FOCUS_0_EOQ', conf_observablesList)
                obs_y_FOCUS_0_EOQ,
            @#end
            @#if ismember('obs_y_FOCUS_20_EOQ', conf_observablesList)
                obs_y_FOCUS_20_EOQ,
            @#end
            @#if ismember('obs_y_FOCUS_25_EOQ', conf_observablesList)
                obs_y_FOCUS_25_EOQ,
            @#end
            @#if ismember('obs_y_FOCUS_50_EOQ', conf_observablesList)
                obs_y_FOCUS_50_EOQ,
            @#end
            @#if ismember('obs_y_FOCUS_75_EOQ', conf_observablesList)
                obs_y_FOCUS_75_EOQ,
            @#end
            @#if ismember('obs_y_FOCUS_80_EOQ', conf_observablesList)
                obs_y_FOCUS_80_EOQ,
            @#end
            @#if ismember('obs_y_FOCUS_100_EOQ', conf_observablesList)
                obs_y_FOCUS_100_EOQ,
            @#end

            @#if ismember('obs_debtOutput_new', conf_observablesList)
                obs_debtOutput_new,
            @#end
             @#if ismember('obs_debtOutput_old', conf_observablesList)
                obs_debtOutput_old,
            @#end
            @#if ismember('obs_tau', conf_observablesList)
                obs_tau,
            @#end
        @#end

        %obs_pii
        %obs_nr

        %obs_pii_FOCUS_median
        %obs_y_FOCUS_median
        %obs_y
        %obs_c obs_y
    @#end

%steady_state_model
%	@#include "closedWithGovernment_steadyState.rs"


@#if conf_estimateModel && ~isempty(conf_observablesList)
    %@#include "paper2_params_estimation_Bayesian.rs"
    %@#include "paper2_params_estimation.rs"
@#end

%%%%%% Specify the planner's objective
%planner_objective{discount = 0.99} -.5*(1*PAI^2+.3*Y^2+0.9*DI^2);
%planner_objective - 1/2 * (1*(Pii-(steady_state(Pii)))^2 + 0.3*(y-steady_state(y))^2);
%planner_objective{discount = beta,commitment=1} - 1/2 * (1*(Pii-(1+piiBar))^2 + 1*(y-ySS)^2);
%planner_objective{discount = beta} -(welfare);
%planner_objective(discount=beta, commitment=1} -(welfare);

%parameterization
% estimate the monetary policy parameter 
%    %paramName, initial, start, end;
%    phi          ,1.5000, 1.0000, 3.0000, normal_pdf(.9);	 
%    gammaTau   	 ,0.1500, 0.0000, 0.2000, normal_pdf(.9);	
