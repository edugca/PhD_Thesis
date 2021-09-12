%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'smoothingCapital';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('smoothingCapital.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eepsM'};
M_.exo_names_tex(1) = {'{\varepsilon^{m}}'};
M_.exo_names_long(1) = {'Monetary shock'};
M_.endo_names = cell(10,1);
M_.endo_names_tex = cell(10,1);
M_.endo_names_long = cell(10,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y_t'};
M_.endo_names_long(1) = {'Output'};
M_.endo_names(2) = {'c'};
M_.endo_names_tex(2) = {'c'};
M_.endo_names_long(2) = {'c'};
M_.endo_names(3) = {'kk'};
M_.endo_names_tex(3) = {'kk'};
M_.endo_names_long(3) = {'kk'};
M_.endo_names(4) = {'kk_lag'};
M_.endo_names_tex(4) = {'kk\_lag'};
M_.endo_names_long(4) = {'kk_lag'};
M_.endo_names(5) = {'nr'};
M_.endo_names_tex(5) = {'i_t'};
M_.endo_names_long(5) = {'Nominal policy rate'};
M_.endo_names(6) = {'r'};
M_.endo_names_tex(6) = {'r_t'};
M_.endo_names_long(6) = {'Real interest rate'};
M_.endo_names(7) = {'ppi'};
M_.endo_names_tex(7) = {'\pi_t'};
M_.endo_names_long(7) = {'Inflation'};
M_.endo_names(8) = {'g'};
M_.endo_names_tex(8) = {'g_t'};
M_.endo_names_long(8) = {'Capital gain'};
M_.endo_names(9) = {'xxiM'};
M_.endo_names_tex(9) = {'\xi^m_t'};
M_.endo_names_long(9) = {'Monetary shock process'};
M_.endo_names(10) = {'AUX_ENDO_LAG_2_1'};
M_.endo_names_tex(10) = {'AUX\_ENDO\_LAG\_2\_1'};
M_.endo_names_long(10) = {'AUX_ENDO_LAG_2_1'};
M_.endo_partitions = struct();
M_.param_names = cell(15,1);
M_.param_names_tex = cell(15,1);
M_.param_names_long = cell(15,1);
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
M_.param_names_long(6) = {'feedback Taylor rule inflation'};
M_.param_names(7) = {'rrhoNr'};
M_.param_names_tex(7) = {'\nu^i'};
M_.param_names_long(7) = {'interest-rate-smoothing'};
M_.param_names(8) = {'rrhoM'};
M_.param_names_tex(8) = {'\rho_d'};
M_.param_names_long(8) = {'autocorrelation monetary'};
M_.param_names(9) = {'ssigmaM'};
M_.param_names_tex(9) = {'\rho_d'};
M_.param_names_long(9) = {'standard deviation monetary'};
M_.param_names(10) = {'ddelta'};
M_.param_names_tex(10) = {'ddelta'};
M_.param_names_long(10) = {'ddelta'};
M_.param_names(11) = {'aalpha'};
M_.param_names_tex(11) = {'aalpha'};
M_.param_names_long(11) = {'aalpha'};
M_.param_names(12) = {'kkappa'};
M_.param_names_tex(12) = {'kkappa'};
M_.param_names_long(12) = {'kkappa'};
M_.param_names(13) = {'ySS'};
M_.param_names_tex(13) = {'ySS'};
M_.param_names_long(13) = {'ySS'};
M_.param_names(14) = {'kkSS'};
M_.param_names_tex(14) = {'kkSS'};
M_.param_names_long(14) = {'kkSS'};
M_.param_names(15) = {'qSS'};
M_.param_names_tex(15) = {'qSS'};
M_.param_names_long(15) = {'qSS'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 10;
M_.param_nbr = 15;
M_.orig_endo_nbr = 9;
M_.aux_vars(1).endo_index = 10;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 3;
M_.aux_vars(1).orig_lead_lag = -1;
M_.aux_vars(1).orig_expr = 'kk(-1)';
M_.predetermined_variables = [ 3 ];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
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
M_.orig_eq_nbr = 9;
M_.eq_nbr = 10;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 2;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 2;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 2;
M_.lead_lag_incidence = [
 0 5 15;
 0 6 16;
 1 7 0;
 0 8 0;
 2 9 0;
 0 10 0;
 0 11 17;
 0 12 18;
 3 13 0;
 4 14 0;]';
M_.nstatic = 2;
M_.nfwrd   = 4;
M_.npred   = 4;
M_.nboth   = 0;
M_.nsfwrd   = 4;
M_.nspred   = 4;
M_.ndynamic   = 8;
M_.dynamic_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'ppi' ;
  4 , 'name' , 'y' ;
  5 , 'name' , 'g' ;
  6 , 'name' , 'nr' ;
  7 , 'name' , '7' ;
  8 , 'name' , 'kk_lag' ;
  9 , 'name' , 'xxiM' ;
};
M_.mapping.y.eqidx = [2 3 4 ];
M_.mapping.c.eqidx = [1 2 3 4 ];
M_.mapping.kk.eqidx = [2 3 4 5 8 ];
M_.mapping.kk_lag.eqidx = [8 ];
M_.mapping.nr.eqidx = [1 6 7 ];
M_.mapping.r.eqidx = [7 ];
M_.mapping.ppi.eqidx = [1 3 6 7 ];
M_.mapping.g.eqidx = [2 5 ];
M_.mapping.xxiM.eqidx = [1 6 9 ];
M_.mapping.eepsM.eqidx = [9 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 5 9 10 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(15, 1);
M_.endo_trends = struct('deflator', cell(10, 1), 'log_deflator', cell(10, 1), 'growth_factor', cell(10, 1), 'log_growth_factor', cell(10, 1));
M_.NNZDerivatives = [38; 0; -1; ];
M_.static_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(2) = 0.99;
bbeta = M_.params(2);
M_.params(3) = 1.00;
eeta = M_.params(3);
M_.params(5) = 0.70;
ttheta = M_.params(5);
M_.params(6) = 1.50;
nnu = M_.params(6);
M_.params(10) = 0.025;
ddelta = M_.params(10);
M_.params(11) = 0.30;
aalpha = M_.params(11);
M_.params(1) = 1.00;
ssigma = M_.params(1);
M_.params(4) = 0.83;
eepsilon = M_.params(4);
M_.params(8) = 0.0;
rrhoM = M_.params(8);
M_.params(7) = 0.0;
rrhoNr = M_.params(7);
M_.params(9) = 1;
ssigmaM = M_.params(9);
M_.params(13) = 1;
ySS = M_.params(13);
M_.params(14) = 5.5;
kkSS = M_.params(14);
M_.params(15) = 1;
qSS = M_.params(15);
M_.params(12) = 0;
kkappa = M_.params(12);
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.01)^2;
write_latex_parameter_table;
write_latex_definitions;
options_.TeX = true;
options_.impulse_responses.plot_threshold = 1e-10;
options_.irf = 20;
options_.nograph = true;
options_.order = 1;
var_list_ = {'xxiM';'kk_lag';'y';'c';'r';'ppi';'nr'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
collect_latex_files;
save('smoothingCapital_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('smoothingCapital_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('smoothingCapital_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('smoothingCapital_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('smoothingCapital_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('smoothingCapital_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('smoothingCapital_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
disp('Note: 1 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
