%%%%%%%%% Amaral(2021): New-Keynesian Models with Defaultables Policy Assets in Zero-Net-Supply

@#define flexiblePrices = false

var
    y   $y_t$ (long_name='Output')
    nr  $i_t$ (long_name='Nominal policy rate')
    r   $r_t$ (long_name='Real interest rate')
    ppi $\pi_t$ (long_name='Inflation')
    iiota $y_t$ (long_name='Policy rule intercept')


    %yRF   $y^{RF}_t$ (long_name='RF Output')
    nrRF  $i^{RF}_t$ (long_name='RF Nominal policy rate')
    rRF   $r^{RF}_t$ (long_name='RF Real interest rate')
    %ppiRF $\pi^{RF}_t$ (long_name='RF Inflation')

    nrPremium  $i^{Premium}_t$ (long_name='Default premium')

    defExp $\mathcal{D}_t$ (long_name='Expected default probability')
    recExp $\omega_t$ (long_name='Expected recovery rate')

    xxiM $\xi^m_t$ (long_name='Monetary shock process')
    xxiD $\xi^m_t$ (long_name='Default shock process')
    xxiR $\xi^m_t$ (long_name='Recovery rate shock process')

    rDirEff $\xi^m_t$ (long_name='Direct effect of monentary shock to the real rate')
    rIndEff $\xi^m_t$ (long_name='Indirect effect of monentary shock to the real rate')
    rExpEff $\xi^m_t$ (long_name='Expectations effect of monentary shock to the real rate')

    %obs_y
    %obs_pii
    %obs_swap_PreDI_3m
    %obs_defLC
    ;

varexo
    eepsM ${\varepsilon^{m}}$ (long_name='Monetary shock')
    eepsD ${\varepsilon^{d}}$ (long_name='Default shock')
    eepsR ${\varepsilon^{r}}$ (long_name='Recovery rate shock')
    
    %eepsPPI
    %eepsY
;

parameters
    ssigma      $\sigma$ (long_name='risk aversion')
    bbeta       $\beta$ (long_name='discount factor')
    eeta        $\eta$ (long_name='Inverse Frisch elasticity')
    eepsilon    $\varepsilon$ (long_name='Intermediate goods elasticity')
    ttheta      $\theta$ (long_name='Calvo parameter')    
    nnu         $\nu$ (long_name='Feedback Taylor rule inflation')
    
    rrhoM       $\rho_d$ (long_name='autocorrelation monetary')
    rrhoD       $\rho_d$ (long_name='autocorrelation default')
    rrhoR       $\rho_d$ (long_name='autocorrelation recovery')

    ssigmaM       $\rho_d$ (long_name='standard deviation monetary')
    ssigmaD       $\rho_d$ (long_name='standard deviation default')
    ssigmaR       $\rho_d$ (long_name='standard deviation recovery')

    defSS       $\rho_d$ (long_name='steady-state default probability')
    recSS       $\rho_d$ (long_name='steady-state recovery rate')
    nrSS

    ytrend
    ppitrend
    nrtrend
    nnuY         $\nu_y$ (long_name='Feedback Taylor rule output gap')

    nnuNr
;

% set deep parameters
eeta        = 5.00;
ttheta      = 0.75; %0.75
bbeta       = 0.99;
defSS       = 0.01;
recSS       = 0.40;
eepsilon    = 6.00;
nnu         = 1.01;
nnuNr       = 0;
rrhoM       = 0.50;
rrhoD       = 0.50;
rrhoR       = 0.50;
ssigmaM     = 0.01;
ssigmaD     = 0.01;
ssigmaR     = 0.01;

piiSS = 0;
rSS = 1/bbeta - 1;
iiotaSS = -1 + ( (1+rSS)*(1+piiSS)-defSS*recSS )/(1-defSS) ;
nrSS = iiotaSS;

// Not being used
ssigma  = 1;
vvarphi = 5;
aalpha  = 1/4;
upsilon = 0.4;

// Estimation
ytrend      = 0;
ppitrend    = 0;
nrtrend     = 0;
nnuY        = 0;

model(linear);

# oomega = (1+eeta)*(1-ttheta)*(1-ttheta*bbeta) / ttheta ;
# pphii = 1 - bbeta*defSS*recSS ;
# pphid = (-defSS + bbeta*defSS*recSS) / (1 - defSS) ;
# pphir = bbeta * defSS * recSS ;

# cM = -pphii*(1-bbeta*rrhoM) / (oomega*(pphii*nnu - rrhoM) + (1-bbeta*rrhoM)*(1-rrhoM) ) ;
# cD = -pphid*(1-bbeta*rrhoD) / (oomega*(pphii*nnu - rrhoD) + (1-bbeta*rrhoD)*(1-rrhoD) ) ;
# cR = -pphir*(1-bbeta*rrhoR) / (oomega*(pphii*nnu - rrhoR) + (1-bbeta*rrhoR)*(1-rrhoR) ) ;

# dM = -pphii*oomega / (oomega*(pphii*nnu - rrhoM) + (1-bbeta*rrhoM)*(1-rrhoM) ) ;
# dD = -pphid*oomega / (oomega*(pphii*nnu - rrhoD) + (1-bbeta*rrhoD)*(1-rrhoD) ) ;
# dR = -pphir*oomega / (oomega*(pphii*nnu - rrhoR) + (1-bbeta*rrhoR)*(1-rrhoR) ) ;

# ggammaI = (1-defSS)*nrSS/(1+nrSS) ;
# ggammaD = (-1-nrSS+recSS)*defSS/(1+nrSS) ;
# ggammaR = (defSS*recSS)/(1+nrSS) ;

@#if flexiblePrices
    y = 0 ;
@#else
    ppi = oomega*y + bbeta*ppi(+1) ;
@#endif

y = y(+1) - (pphii*nr - ppi(+1) + pphid*defExp + pphir*recExp) ;
y = y(+1) - (nrRF - ppi(+1)) ;

nr = nnuNr*nr(-1) + (1-nnuNr)* (iiota + nnu*ppi + nnuY*y) + xxiM ;
iiota = 0 ;

%r = rRF ;
r = -ppi(+1) + ggammaI*nr + ggammaD*defExp + ggammaR*recExp ;


// Risk-free
%yRF = yRF(+1) - (nrRF - ppiRF(+1));
%ppiRF = oomega*yRF + bbeta*ppiRF(+1);
%nrRF = nnu*ppiRF + xxiM;
rRF = nrRF - ppi(+1);

nrPremium = nr - nrRF;

defExp = xxiD;
recExp = xxiR;

// Shocks
xxiM = rrhoM*xxiM(-1) + ssigmaM * eepsM;
xxiD = rrhoD*xxiD(-1) + ssigmaD * eepsD;
xxiR = rrhoR*xxiR(-1) + ssigmaR * eepsR;

// Decomposition
rDirEff = 1 * xxiM;
rIndEff = nnu * dM * xxiM;
rExpEff = - rrhoM * dM * xxiM;

// Measurement equations
%obs_y                   = y     + ytrend + eepsY ;
%obs_pii                 = ppi   + ppitrend;
%obs_swap_PreDI_3m/100   = nr    + nrtrend;
%obs_defLC/4             = defExp   + defSS; 

end;

shocks;
    var eepsM           = ssigmaM^2; //standard deviation of the monetary shock
    var eepsD           = ssigmaD^2; //standard deviation of the default shock
    var eepsR           = ssigmaR^2; //standard deviation of the recovery rate shock
    
    %var eepsPPI         = 0.01; //standard deviation of inflation
    %var eepsY         = 0.01; //standard deviation of inflation

    %var eps_a, eps_star = 0.3*0.0071*0.0078; //covariance, constructed from correlation and standard deviations
end;

//generate LaTeX-files with equations and parameterization
write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

% nograph
% irf_plot_threshold=1e-10
stoch_simul(TeX,order=1,irf=20) xxiM xxiD xxiR y ppi nr r rRF defExp recExp;

collect_latex_files;

