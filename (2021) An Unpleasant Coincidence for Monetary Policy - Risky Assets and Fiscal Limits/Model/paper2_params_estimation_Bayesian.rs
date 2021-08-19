% Parameters that enter into the calculation of fiscal limits
%aTildeSS, recTildeSS, gSS, defSS, ySS, zSS, tLSSS, yMaxSS, tauMaxSS, ...
%beta, chi, eta, kkSS, sigma, alphaG, gammaGPSI, fracNR, ...
%rhoA, rhoR, rhoD, rhoGG, rhoGY, ...
%sigmaA, sigmaR, sigmaD, sigmaG, ...


%%%%%% RESTRICTIONS ON ESTIMATION
parameter_restrictions;

% Parameters that enter into the calculation of fiscal limits
%aTildeSS, recTildeSS, gSS, defSS, ySS, zSS, tLSSS, yMaxSS, tauMaxSS, ...
%beta, chi, eta, kkSS, sigma, alphaG, gammaGPSI, fracNR, ...
%rhoA, rhoR, rhoD, rhoGG, rhoGY, ...
%sigmaA, sigmaR, sigmaD, sigmaG, ...

@#if conf_prodProcessType == "Dual_LogNormal"
    %defSS <= 0.5;
    aTildeSS > 0;
    recTildeSS > 0;
    aTildeSS > recTildeSS + 0.01;
@#end

%%%%%% PARAMETERS FOR ESTIMATION
% name, starting value, lower quantile, upper quantile, dist(% between quantiles)
parameterization;

    @#if conf_occBinConstraint == "FiscalLimit"
        %fisLim_tp_1_2,  0.5,    0.0001, 0.1, beta_pdf(.9);
        %fisLim_tp_2_1,  1.0,    0.9, 0.9999, beta_pdf(.9);
        %taxLim_tp_1_2,  0.5,    0.0001, 0.1, beta_pdf(.9);
        %taxLim_tp_2_1,  1.0,    0.9, 0.9999, beta_pdf(.9);
    @#end

    @#if conf_defRProcessType == "Null"
        defSS,  0.0;
        rhoD,   0.0;
        sigmaD, 0.0;
    @#elseif conf_defRProcessType == "Fixed"
        defSS,  0.5;
        %defSS,  0.25,    0.001, 0.5, normal_pdf(1.00);
        rhoD,   0.0;
        sigmaD, 0.0;
    @#elseif conf_defRProcessType == "LogNormal"
        %defSS, 0.5;
        defSS,  0.5,    0.2, 0.5, beta_pdf(.9);
        rhoD,   0.5,    0.4, 0.8, beta_pdf(.9);
        sigmaD, 0.01,   0.00001, 0.3, inv_gamma_pdf(0.999);
    @#end

    @#if conf_prodProcessType == "MarkovSwitching"
        sigmaD(s,1), 0.005,   0.0005, 0.1, inv_gamma_pdf(0.999);
        sigmaD(s,2), 0.005,   0.0005, 0.1, inv_gamma_pdf(0.999);
        s_tp_1_2, 1-0.75;
        s_tp_2_1, 1-0.5;
    @#end
    
    @#if (conf_shockOnEffectiveProductivity == "OnlyExpectation") | (conf_shockOnEffectiveProductivity == "OnlyRealization")
        sigmaPsi, 0.01,   0.0001, 0.1, inv_gamma_pdf(0.999);
    @#end

    @#if conf_kVariable == "Fixed"
        sigmaK, 0;
    @#else
        sigmaK, 0.008,   0.01, 0.1, inv_gamma_pdf(.99);
    @#end

    %%% PREFERENCES
    %beta,   0.99;
    sigma,  1.5,     1.5, 3.5,    gamma_pdf(.99);
    %eta,    0.25,   0.1, 5.0, inv_gamma_pdf(.99);
    %chi,    1,     1.0, 2.0,     gamma_pdf(0.99); %%%FISCAL LIMITS ARE VERY SENSITIVE TO IT
    
    
    rhoA,   0.5,    0.25,    0.75, beta_pdf(.9);
    sigmaA, 0.01,    0.0001, 0.035, inv_gamma_pdf(0.999);
    %ytrend,    0.005,  0.005,    0.01,  normal_pdf(.9);

    @#if ismember('obs_c', conf_observablesList)
        rhoBeta,   0.5,    0.25,    0.75, beta_pdf(.9);
        sigmaBeta, 0.01,    0.0001, 0.02, inv_gamma_pdf(0.999);
    @#end

    @#if (sum(ismember(conf_observablesList, 'obs_pii')) > 0) | ismember('obs_nr', conf_observablesList)
        
        %rhoM, 0.5;
        rhoM,   0.5,    0.25,    0.75, beta_pdf(.9);
        sigmaM, 0.01,    0.0001, 0.02, inv_gamma_pdf(0.999);

        %piiBar, 0.01,   0,    0.02,  normal_pdf(.95);
        phiNr,    0.75,  0.5,    0.99,    beta_pdf(.95);
        %phi_Exp,  3,   1.001,  3.5,    normal_pdf(.99);
        phi_Y,    0.0001,      0.0,    1.0,    normal_pdf(.95);
        %phi_dY,   0.0001,      0.0,    1.0,    normal_pdf(.95);

        @#if conf_occBinConstraint == "None"
            phi,    3,   2,    3.5,  normal_pdf(.95);
        @#elseif conf_occBinConstraint == "FiscalLimit"
            phi,    3,   2,    3.5,  normal_pdf(.95);
            %phi(fisLim,1),    3,   2,    3.5,  normal_pdf(.95);
        @#end
    @#end

    @#if sum(ismember(conf_observablesList, 'obs_defLC')) > 0
        rhoPolDef,      0.15,    0.25,    0.75, beta_pdf(.9);
        sigmaPolDef,    0.01,    0.0001, 0.02, inv_gamma_pdf(0.999);
        polDefWedge,    0.01,   0.0001,    0.1,  normal_pdf(.95);

        %muFisLim_Intercept,    0.0,   -5, 5,  normal_pdf(.9);
        %muFisLim_aTilde,       0.0,   -5, 5,  normal_pdf(.9);
        %muFisLim_g,            0.0,   -5, 5,  normal_pdf(.9);
        %muFisLim_bbeta,        0.0,   -5, 5,  normal_pdf(.9);

        %stdFisLim_Intercept,    0.0,   -5, 5,  normal_pdf(.9);
        %stdFisLim_aTilde,       0.0,   -5, 5,  normal_pdf(.9);
        %stdFisLim_g,            0.0,   -5, 5,  normal_pdf(.9);
        %stdFisLim_bbeta,        0.0,   -5, 5,  normal_pdf(.9);

        %probDefFisLim_Param_1,       0.0,   -5, 5,  normal_pdf(.9);
        %probDefFisLim_Param_2,       0.0,   -5, 5,  normal_pdf(.9);
        %probDefFisLim_Param_3,       0.0,   -5, 5,  normal_pdf(.9);
        %probDefFisLim_Param_B,       0.0,   -5, 5,  normal_pdf(.9);
        %probDefFisLim_Param_A,       0.0,   -5, 5,  normal_pdf(.9);
        %probDefFisLim_Param_G,       0.0,   -5, 5,  normal_pdf(.9);
        %probDefFisLim_Param_Beta,    0.0,   -5, 5,  normal_pdf(.9);
    @#end
    
    @#if conf_prodProcessType == "Single_LogNormal"
        rhoR,   0.0;
        sigmaR, 0.00;
    @#elseif conf_prodProcessType == "Dual_LogNormal"
        aTildeSS,  0.20,  0.0001,    0.25,  normal_pdf(0.75);

        rhoR,   0.9,    0.25,    0.75, beta_pdf(.9);
        sigmaR, 0.01,    0.0001, 0.025, inv_gamma_pdf(0.999);
        recTildeSS,  0.15,  0.0001,    0.25,  normal_pdf(0.75);
    @#end

    %kBar,   1;
    %pBar,   1;

    %%%%%%%%%%%%%%%% Labor
    @#if ismember('obs_unempRate', conf_observablesList)
        unempSS,    0.045,   0.01, 0.06, normal_pdf(0.95);
    @#end

    %%%%%%%%%%%%%%%% Government
    @#if ismember('obs_g', conf_observablesList)
        rhoGG,      0.7,    0.25,    0.75, beta_pdf(.9);
        rhoGY,      0.1,    0.25,    0.75, beta_pdf(.9);
        %rhoTau,    0.7,    0.3, 0.7, beta_pdf(.9);
        %sigmaTau,     0.0045, 0.0212, 0.2810, inv_gamma_pdf(0.999);
        sigmaG,     0.01,    0.0001, 0.035, inv_gamma_pdf(0.999);
        alphaG,     0,      -1, 1,      normal_pdf(0.95);
        gammaGPSI,  0,      -0.25, 0.5,  normal_pdf(0.95);

        %gammaTau(taxLim,1),    2.0,   1.0, 3.0, normal_pdf(0.9);
    @#end

    %%%%%%%%%%%%%%%% Turn Off shocks
    @#if ismember('sigmaTau', conf_shocksTurnedOff)
        sigmaTau, 0;
    @#end
    %%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%% MEASUREMENT ERRORS

    @#if conf_measurementErrors
        @#if ismember('obs_y', conf_measurementErrorsList)
            stderr_obs_y,   0.01,    0.0001, 0.04, inv_gamma_pdf(0.999);
        @#end

         @#if ismember('obs_defLC', conf_measurementErrorsList)
            stderr_obs_defLC,   0.01,    0.0001, 0.04, inv_gamma_pdf(0.999);
        @#end

        @#if ismember('obs_unempRate', conf_measurementErrorsList)
            stderr_obs_unempRate,   0.01,   0.0001, 0.1, inv_gamma_pdf(0.999);
        @#end

        @#if ismember('obs_wage', conf_measurementErrorsList)
            stderr_obs_wage,   0.01,    0.0001, 0.04, inv_gamma_pdf(0.999);
        @#end

        @#if ismember('obs_pii', conf_measurementErrorsList)
            stderr_obs_pii,   0.01,   0.0001, 0.1, inv_gamma_pdf(0.999);
        @#end

        @#if ismember('obs_nr', conf_measurementErrorsList)
            stderr_obs_nr,   0.0025,   0.001, 0.05, inv_gamma_pdf(.99);
        @#end

        @#if ismember('obs_swap_PreDI_3m', conf_measurementErrorsList)
            stderr_obs_swap_PreDI_3m,   0.0025,   0.001, 0.05, inv_gamma_pdf(.99);
        @#end
        
        @#if ismember('obs_pii_FOCUS_Median', conf_measurementErrorsList)
            stderr_obs_pii_FOCUS_Median,   0.005,   0.0001, 0.03, inv_gamma_pdf(0.999);
        @#end

        @#if ismember('obs_recTilde_FOCUS_0', conf_measurementErrorsList)
            stderr_obs_recTilde_FOCUS_0,   0.001,   0.001, 0.05, inv_gamma_pdf(.99);
        @#end

        @#if ismember('obs_recTilde_FOCUS_25', conf_measurementErrorsList)
            stderr_obs_recTilde_FOCUS_25,   0.001,   0.001, 0.05, inv_gamma_pdf(.99);
        @#end

        @#if ismember('obs_aTilde_FOCUS_75', conf_measurementErrorsList)
            stderr_obs_aTilde_FOCUS_75,   0.001,   0.001, 0.05, inv_gamma_pdf(.99);
        @#end

        @#if ismember('obs_aTilde_FOCUS_100', conf_measurementErrorsList)
            stderr_obs_aTilde_FOCUS_100,   0.001,   0.001, 0.05, inv_gamma_pdf(.99);
        @#end

        @#if ismember('obs_y_FOCUS_20_EOQ', conf_measurementErrorsList)
            stderr_obs_y_FOCUS_20_EOQ,   0.001,   0.001, 0.05, inv_gamma_pdf(.95);
        @#end

        @#if ismember('obs_y_FOCUS_80_EOQ', conf_measurementErrorsList)
            stderr_obs_y_FOCUS_80_EOQ,   0.001,   0.001, 0.05, inv_gamma_pdf(.95);
        @#end

    @#end

    %stderr_obs_yPotentialQoQ, 0.01,   0.01, 0.1, inv_gamma_pdf(.99);
    %stderr_obs_lPotentialQoQ, 0.01,   0.01, 0.1, inv_gamma_pdf(.99);
    @#if conf_kVariable == "Exogenous"
        stderr_obs_kPotentialQoQ, 0.01,   0.001, 0.05, inv_gamma_pdf(.99);
    @#end