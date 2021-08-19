% Parameters that enter into the calculation of fiscal limits
%aTildeSS, recTildeSS, gSS, defSS, ySS, zSS, tLSSS, yMaxSS, tauMaxSS, ...
%beta, chi, eta, kkSS, sigma, alphaG, gammaGPSI, fracNR, ...
%rhoA, rhoR, rhoD, rhoGG, rhoGY, ...
%sigmaA, sigmaR, sigmaD, sigmaG, ...


%%%%%% RESTRICTIONS ON ESTIMATION
parameter_restrictions;
@#if conf_prodProcessType == "Dual_LogNormal"
    %defSS <= 0.5;
    aTildeSS > 0;
    recTildeSS > 0;
    aTildeSS > recTildeSS + 0.01;
@#end

%@#if conf_occBinConstraint == "None"
%    phi_fisLim_2 = phi_fisLim_1
%@#end

%%%%%% PARAMETERS FOR ESTIMATION
parameterization;

%%% PREFERENCES
%beta,   0.989, 0.95,    0.999;
%chi,    1,      0.1,     5;
sigma,    1.5,      1,   5;

rhoBeta,   0.5,    0.1, 0.9;
sigmaBeta, 0.01,   0.00001, 0.03;

%%% PRODUCTIVITY
rhoA,   0.884,    0.3, 0.99;
sigmaA, 0.0125,   0.00005, 0.05;
%ytrend,    0.005,  -0.015,    0.015;

%%% INFLATION AND PREFERENCES
@#if (sum(ismember(conf_observablesList, 'obs_pii')) > 0) || ismember('obs_nr', conf_observablesList)
    rhoM,   0.5,    0.1, 0.7;
    sigmaM, 0.005,   0.00001, 0.03;
    %piiBar,     0.011,   0,    0.02;

    phiNr,      0.79,  0.5,    0.9;
    %phi_Exp,   3,   1.001,  3.5;
    phi_Y,      0.0001,      0,    1.0;
    %phi_dY,     0.0001,      0,    1.0;

    @#if conf_occBinConstraint == "None"
        phi,    2.43,   1.1,    5.0;
    @#elseif conf_occBinConstraint == "FiscalLimit"
        phi(fisLim,1),    3,   1.001,    3.5;
        %phi(fisLim,2),    3,   1.001,    3.5;
    @#end
@#end

%%% Government
@#if ismember('obs_g', conf_observablesList)
    rhoGG,      0.986,      0.1, 0.99;
    rhoGY,      0.786,      0.1, 0.9;
    
    sigmaG,     0.0045, 0.00001, 0.05;
    alphaG,     0.5,      0, 3;
    gammaGPSI,  0.2,    -1.5, 1.5;

    %rhoPolDef,   0;
    %sigmaPolDef, 0;
    %rhoTau,    0.7,        0.3, 0.7, beta_pdf(.9);
    %gammaTau(taxLim,1),    2.0,   1.0, 3.0, normal_pdf(0.9);
@#end

%%%%%%%%%%%%%%%% MEASUREMENT ERRORS
@#if conf_measurementErrors
    @#if ismember('obs_y', conf_measurementErrorsList)
        stderr_obs_y,   0.01,   0.0001, 0.1;
    @#end

    @#if ismember('obs_unempRate', conf_measurementErrorsList)
        stderr_obs_unempRate,   0.01,   0.0001, 0.1;
    @#end

    @#if ismember('obs_w', conf_measurementErrorsList)
        stderr_obs_w,   0.01,    0.0001, 0.025, inv_gamma_pdf(0.999);
    @#end

    @#if ismember('obs_pii', conf_measurementErrorsList)
        stderr_obs_pii,   0.01,   0.0001, 0.1;
    @#end

    @#if ismember('obs_nr', conf_measurementErrorsList)
        stderr_obs_nr,   0.0025,   0.001, 0.05, inv_gamma_pdf(.99);
    @#end

    @#if ismember('obs_swap_PreDI_3m', conf_measurementErrorsList)
        stderr_obs_swap_PreDI_3m,   0.0025,   0.001, 0.05, inv_gamma_pdf(.99);
    @#end

    @#if ismember('obs_pii_FOCUS_Median', conf_measurementErrorsList)
        stderr_obs_pii_FOCUS_Median,   0.005,   0.0001, 0.03;
    @#end

    @#if ismember('obs_recTilde_FOCUS_0', conf_measurementErrorsList)
        stderr_obs_recTilde_FOCUS_0,   0.001,   0.001, 0.05;
    @#end

    @#if ismember('obs_recTilde_FOCUS_25', conf_measurementErrorsList)
        stderr_obs_recTilde_FOCUS_25,   0.001,   0.001, 0.05;
    @#end

    @#if ismember('obs_aTilde_FOCUS_75', conf_measurementErrorsList)
        stderr_obs_aTilde_FOCUS_75,   0.001,   0.001, 0.05;
    @#end

    @#if ismember('obs_aTilde_FOCUS_100', conf_measurementErrorsList)
        stderr_obs_aTilde_FOCUS_100,   0.001,   0.001, 0.05;
    @#end

    @#if ismember('obs_y_FOCUS_20_EOQ', conf_measurementErrorsList)
        stderr_obs_y_FOCUS_20_EOQ,   0.001,   0.001, 0.05;
    @#end

    @#if ismember('obs_y_FOCUS_80_EOQ', conf_measurementErrorsList)
        stderr_obs_y_FOCUS_80_EOQ,   0.001,   0.001, 0.05;
    @#end

@#end

%stderr_obs_yPotentialQoQ, 0.01,   0.01, 0.1;
%stderr_obs_lPotentialQoQ, 0.01,   0.01, 0.1;
@#if conf_kVariable == "Exogenous"
    stderr_obs_kPotentialQoQ, 0.01,   0.001, 0.05;
@#end