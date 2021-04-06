%%%%%%%%% Amaral(2021): On the Mechanics of New-Keynesian Models: Smoothing the Capital Controversy Out
%%%%%%%%% Replication of Rupert and Sustek (2019) with the addition of Interest-rate-smoothing

var
    y   $y_t$ (long_name='Output')
    c
    kk
    kk_lag
    nr  $i_t$ (long_name='Nominal policy rate')
    r   $r_t$ (long_name='Real interest rate')
    ppi $\pi_t$ (long_name='Inflation')

    g $g_t$  (long_name='Capital gain')
    %q $q_t$  (long_name='Tobin q')

    xxiM $\xi^m_t$ (long_name='Monetary shock process')

    %rDirEff $\xi^m_t$ (long_name='Direct effect of monentary shock to the real rate')
    %rIndEff $\xi^m_t$ (long_name='Indirect effect of monentary shock to the real rate')
    %rExpEff $\xi^m_t$ (long_name='Expectations effect of monentary shock to the real rate')
    ;
predetermined_variables kk;

varexo
    eepsM ${\varepsilon^{m}}$ (long_name='Monetary shock')
;

parameters
    ssigma      $\sigma$ (long_name='risk aversion')
    bbeta       $\beta$ (long_name='discount factor')
    eeta        $\eta$ (long_name='Inverse Frisch elasticity')
    eepsilon    $\varepsilon$ (long_name='Intermediate goods elasticity')
    ttheta      $\theta$ (long_name='Calvo parameter')    
    nnu         $\nu$ (long_name='feedback Taylor rule inflation')

    rrhoNr       $\nu^i$ (long_name='interest-rate-smoothing')

    rrhoM       $\rho_d$ (long_name='autocorrelation monetary')
    ssigmaM     $\rho_d$ (long_name='standard deviation monetary')

    ddelta
    aalpha

    kkappa

    ySS
    kkSS
    qSS
;

% set deep parameters
bbeta       = 0.99;
eeta        = 1.00;
ttheta      = 0.70; % 0.70
nnu         = 1.50;
ddelta      = 0.025; % 0.025
aalpha      = 0.30;
ssigma      = 1.00;
eepsilon    = 0.83;

rrhoM       = 0.0; % 0, 0.1, 0.5, 0.95
rrhoNr      = 0.0; % 0.6
ssigmaM     = 1;

% set steady-state parameters
ySS     = 1;
kkSS    = 5.5; % 5.5 the closest
qSS     = 1;
kkappa  = 0.0; %capital adjustment costs

model(linear);

# kkappaLine    = kkappa*kkSS;
# cSS           = ySS -ddelta*kkSS;
# wSS           = cSS*(ySS/(kkSS^aalpha))^(eeta/(1-aalpha));
# rSS           = 1/bbeta - 1 + ddelta;
# cchiSS        = wSS/(1-aalpha)*(ySS/kkSS)^(aalpha/(1-aalpha));
# ppsi          = cchiSS*(1-ttheta)*(1-ttheta*bbeta)/ttheta;

-c = -c(+1) + rrhoNr*nr(-1) + (1-rrhoNr)* (nnu*ppi) - ppi(+1) + ssigmaM*xxiM ;
-c = -c(+1) + g(+1) + rSS*( c(+1) + (1+eeta)/(1-aalpha)*y(+1) - (1+aalpha*eeta)/(1-aalpha)*kk(+1) ) ;
ppi = ppsi*( (eeta + aalpha)/(1 - aalpha)*y - (aalpha*(1+eeta))/(1-aalpha)*kk + c ) + bbeta*ppi(+1) ;
y = cSS/ySS*c + kkSS/ySS*kk(+1) - (1-ddelta)*kkSS/ySS*kk ;

g = kkappaLine*(kk(+1) - kk) - kkappaLine*(kk - kk(-1));

nr = rrhoNr*nr(-1) + (1-rrhoNr)* ( nnu*ppi ) + ssigmaM*xxiM ;
nr = r + ppi(+1);

kk_lag = kk(-1);

// Shocks
xxiM = rrhoM*xxiM(-1) + ssigmaM * eepsM;

end;

steady_state_model;

ppi = 0;
c   = 0;
y   = 0;
kk = 0;
nr = 0;
r = 0;
q = 0;
g = 0;
xxiM = 0;

end;

steady;

shocks;
    var eepsM           = 1^2; //standard deviation of the monetary shock

    %var eps_a, eps_star = 0.3*0.0071*0.0078; //covariance, constructed from correlation and standard deviations
end;

//generate LaTeX-files with equations and parameterization
write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

% nograph
stoch_simul(TeX,order=1,irf=20, irf_plot_threshold=1e-10, nograph) xxiM kk_lag y c r ppi nr ;

collect_latex_files;

