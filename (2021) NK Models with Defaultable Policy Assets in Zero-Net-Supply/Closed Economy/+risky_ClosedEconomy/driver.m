%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'risky_ClosedEconomy';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('risky_ClosedEconomy.log');
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'eepsM'};
M_.exo_names_tex(1) = {'{\varepsilon^{m}}'};
M_.exo_names_long(1) = {'Monetary shock'};
M_.exo_names(2) = {'eepsD'};
M_.exo_names_tex(2) = {'{\varepsilon^{d}}'};
M_.exo_names_long(2) = {'Default shock'};
M_.exo_names(3) = {'eepsR'};
M_.exo_names_tex(3) = {'{\varepsilon^{r}}'};
M_.exo_names_long(3) = {'Recovery rate shock'};
M_.endo_names = cell(16,1);
M_.endo_names_tex = cell(16,1);
M_.endo_names_long = cell(16,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y_t'};
M_.endo_names_long(1) = {'Output'};
M_.endo_names(2) = {'nr'};
M_.endo_names_tex(2) = {'i_t'};
M_.endo_names_long(2) = {'Nominal policy rate'};
M_.endo_names(3) = {'r'};
M_.endo_names_tex(3) = {'r_t'};
M_.endo_names_long(3) = {'Real interest rate'};
M_.endo_names(4) = {'ppi'};
M_.endo_names_tex(4) = {'\pi_t'};
M_.endo_names_long(4) = {'Inflation'};
M_.endo_names(5) = {'iiota'};
M_.endo_names_tex(5) = {'y_t'};
M_.endo_names_long(5) = {'Policy rule intercept'};
M_.endo_names(6) = {'nrRF'};
M_.endo_names_tex(6) = {'i^{RF}_t'};
M_.endo_names_long(6) = {'RF Nominal policy rate'};
M_.endo_names(7) = {'rRF'};
M_.endo_names_tex(7) = {'r^{RF}_t'};
M_.endo_names_long(7) = {'RF Real interest rate'};
M_.endo_names(8) = {'nrPremium'};
M_.endo_names_tex(8) = {'i^{Premium}_t'};
M_.endo_names_long(8) = {'Default premium'};
M_.endo_names(9) = {'defExp'};
M_.endo_names_tex(9) = {'\mathcal{D}_t'};
M_.endo_names_long(9) = {'Expected default probability'};
M_.endo_names(10) = {'recExp'};
M_.endo_names_tex(10) = {'\omega_t'};
M_.endo_names_long(10) = {'Expected recovery rate'};
M_.endo_names(11) = {'xxiM'};
M_.endo_names_tex(11) = {'\xi^m_t'};
M_.endo_names_long(11) = {'Monetary shock process'};
M_.endo_names(12) = {'xxiD'};
M_.endo_names_tex(12) = {'\xi^m_t'};
M_.endo_names_long(12) = {'Default shock process'};
M_.endo_names(13) = {'xxiR'};
M_.endo_names_tex(13) = {'\xi^m_t'};
M_.endo_names_long(13) = {'Recovery rate shock process'};
M_.endo_names(14) = {'rDirEff'};
M_.endo_names_tex(14) = {'\xi^m_t'};
M_.endo_names_long(14) = {'Direct effect of monentary shock to the real rate'};
M_.endo_names(15) = {'rIndEff'};
M_.endo_names_tex(15) = {'\xi^m_t'};
M_.endo_names_long(15) = {'Indirect effect of monentary shock to the real rate'};
M_.endo_names(16) = {'rExpEff'};
M_.endo_names_tex(16) = {'\xi^m_t'};
M_.endo_names_long(16) = {'Expectations effect of monentary shock to the real rate'};
M_.endo_partitions = struct();
M_.param_names = cell(20,1);
M_.param_names_tex = cell(20,1);
M_.param_names_long = cell(20,1);
M_.param_names(1) = {'ssigma'};
M_.param_names_tex(1) = {'\sigma'};
M_.param_names_long(1) = {'risk aversion'};
M_.param_names(2) = {'bbeta'};
M_.param_names_tex(2) = {'\beta'};
M_.param_names_long(2) = {'discount factor'};
M_.param_names(3) = {'eeta'};
M_.param_names_tex(3) = {'\eta'};
M_.param_names_long(3) = {'Inverse Frisch elasticity'};
M_.param_names(4) = {'eepsilon'};
M_.param_names_tex(4) = {'\varepsilon'};
M_.param_names_long(4) = {'Intermediate goods elasticity'};
M_.param_names(5) = {'ttheta'};
M_.param_names_tex(5) = {'\theta'};
M_.param_names_long(5) = {'Calvo parameter'};
M_.param_names(6) = {'nnu'};
M_.param_names_tex(6) = {'\nu'};
M_.param_names_long(6) = {'Feedback Taylor rule inflation'};
M_.param_names(7) = {'rrhoM'};
M_.param_names_tex(7) = {'\rho_d'};
M_.param_names_long(7) = {'autocorrelation monetary'};
M_.param_names(8) = {'rrhoD'};
M_.param_names_tex(8) = {'\rho_d'};
M_.param_names_long(8) = {'autocorrelation default'};
M_.param_names(9) = {'rrhoR'};
M_.param_names_tex(9) = {'\rho_d'};
M_.param_names_long(9) = {'autocorrelation recovery'};
M_.param_names(10) = {'ssigmaM'};
M_.param_names_tex(10) = {'\rho_d'};
M_.param_names_long(10) = {'standard deviation monetary'};
M_.param_names(11) = {'ssigmaD'};
M_.param_names_tex(11) = {'\rho_d'};
M_.param_names_long(11) = {'standard deviation default'};
M_.param_names(12) = {'ssigmaR'};
M_.param_names_tex(12) = {'\rho_d'};
M_.param_names_long(12) = {'standard deviation recovery'};
M_.param_names(13) = {'defSS'};
M_.param_names_tex(13) = {'\rho_d'};
M_.param_names_long(13) = {'steady-state default probability'};
M_.param_names(14) = {'recSS'};
M_.param_names_tex(14) = {'\rho_d'};
M_.param_names_long(14) = {'steady-state recovery rate'};
M_.param_names(15) = {'nrSS'};
M_.param_names_tex(15) = {'nrSS'};
M_.param_names_long(15) = {'nrSS'};
M_.param_names(16) = {'ytrend'};
M_.param_names_tex(16) = {'ytrend'};
M_.param_names_long(16) = {'ytrend'};
M_.param_names(17) = {'ppitrend'};
M_.param_names_tex(17) = {'ppitrend'};
M_.param_names_long(17) = {'ppitrend'};
M_.param_names(18) = {'nrtrend'};
M_.param_names_tex(18) = {'nrtrend'};
M_.param_names_long(18) = {'nrtrend'};
M_.param_names(19) = {'nnuY'};
M_.param_names_tex(19) = {'\nu_y'};
M_.param_names_long(19) = {'Feedback Taylor rule output gap'};
M_.param_names(20) = {'nnuNr'};
M_.param_names_tex(20) = {'nnuNr'};
M_.param_names_long(20) = {'nnuNr'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 16;
M_.param_nbr = 20;
M_.orig_endo_nbr = 16;
M_.aux_vars = [];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 16;
M_.eq_nbr = 16;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 5 21;
 1 6 0;
 0 7 0;
 0 8 22;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 2 15 0;
 3 16 0;
 4 17 0;
 0 18 0;
 0 19 0;
 0 20 0;]';
M_.nstatic = 10;
M_.nfwrd   = 2;
M_.npred   = 4;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.dynamic_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'ppi' ;
  2 , 'name' , 'y' ;
  3 , 'name' , '3' ;
  4 , 'name' , 'nr' ;
  5 , 'name' , 'iiota' ;
  6 , 'name' , 'r' ;
  7 , 'name' , 'rRF' ;
  8 , 'name' , 'nrPremium' ;
  9 , 'name' , 'defExp' ;
  10 , 'name' , 'recExp' ;
  11 , 'name' , 'xxiM' ;
  12 , 'name' , 'xxiD' ;
  13 , 'name' , 'xxiR' ;
  14 , 'name' , 'rDirEff' ;
  15 , 'name' , 'rIndEff' ;
  16 , 'name' , 'rExpEff' ;
};
M_.mapping.y.eqidx = [1 2 3 4 ];
M_.mapping.nr.eqidx = [2 4 6 8 ];
M_.mapping.r.eqidx = [6 ];
M_.mapping.ppi.eqidx = [1 2 3 4 6 7 ];
M_.mapping.iiota.eqidx = [5 ];
M_.mapping.nrRF.eqidx = [3 7 8 ];
M_.mapping.rRF.eqidx = [7 ];
M_.mapping.nrPremium.eqidx = [8 ];
M_.mapping.defExp.eqidx = [2 6 9 ];
M_.mapping.recExp.eqidx = [2 6 10 ];
M_.mapping.xxiM.eqidx = [4 11 14 15 16 ];
M_.mapping.xxiD.eqidx = [9 12 ];
M_.mapping.xxiR.eqidx = [10 13 ];
M_.mapping.rDirEff.eqidx = [14 ];
M_.mapping.rIndEff.eqidx = [15 ];
M_.mapping.rExpEff.eqidx = [16 ];
M_.mapping.eepsM.eqidx = [11 ];
M_.mapping.eepsD.eqidx = [12 ];
M_.mapping.eepsR.eqidx = [13 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [2 11 12 13 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(16, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(20, 1);
M_.endo_trends = struct('deflator', cell(16, 1), 'log_deflator', cell(16, 1), 'growth_factor', cell(16, 1), 'log_growth_factor', cell(16, 1));
M_.NNZDerivatives = [49; 0; -1; ];
M_.static_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(3) = 5.00;
eeta = M_.params(3);
M_.params(5) = 0.75;
ttheta = M_.params(5);
M_.params(2) = 0.99;
bbeta = M_.params(2);
M_.params(13) = 0.01;
defSS = M_.params(13);
M_.params(14) = 0.40;
recSS = M_.params(14);
M_.params(4) = 6.00;
eepsilon = M_.params(4);
M_.params(6) = 1.01;
nnu = M_.params(6);
M_.params(20) = 0;
nnuNr = M_.params(20);
M_.params(7) = 0.50;
rrhoM = M_.params(7);
M_.params(8) = 0.50;
rrhoD = M_.params(8);
M_.params(9) = 0.50;
rrhoR = M_.params(9);
M_.params(10) = 0.01;
ssigmaM = M_.params(10);
M_.params(11) = 0.01;
ssigmaD = M_.params(11);
M_.params(12) = 0.01;
ssigmaR = M_.params(12);
piiSS = 0;
rSS = 1/bbeta - 1;
iiotaSS = -1 + ( (1+rSS)*(1+piiSS)-defSS*recSS )/(1-defSS) ;
M_.params(15) = iiotaSS;
nrSS = M_.params(15);
M_.params(1) = 1;
ssigma = M_.params(1);
vvarphi = 5;
aalpha  = 1/4;
upsilon = 0.4;
M_.params(16) = 0;
ytrend = M_.params(16);
M_.params(17) = 0;
ppitrend = M_.params(17);
M_.params(18) = 0;
nrtrend = M_.params(18);
M_.params(19) = 0;
nnuY = M_.params(19);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(10)^2;
M_.Sigma_e(2, 2) = M_.params(11)^2;
M_.Sigma_e(3, 3) = M_.params(12)^2;
write_latex_parameter_table;
write_latex_definitions;
options_.TeX = true;
options_.irf = 20;
options_.order = 1;
var_list_ = {'xxiM';'xxiD';'xxiR';'y';'ppi';'nr';'r';'rRF';'defExp';'recExp'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
collect_latex_files;
save('risky_ClosedEconomy_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('risky_ClosedEconomy_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('risky_ClosedEconomy_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('risky_ClosedEconomy_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('risky_ClosedEconomy_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('risky_ClosedEconomy_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('risky_ClosedEconomy_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
