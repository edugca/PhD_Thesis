%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Brault_Khan_JMCB_2019';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Brault_Khan_JMCB_2019.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps_m'};
M_.exo_names_tex(1) = {'eps\_m'};
M_.exo_names_long(1) = {'eps_m'};
M_.endo_names = cell(15,1);
M_.endo_names_tex = cell(15,1);
M_.endo_names_long = cell(15,1);
M_.endo_names(1) = {'Y'};
M_.endo_names_tex(1) = {'Y'};
M_.endo_names_long(1) = {'Y'};
M_.endo_names(2) = {'C'};
M_.endo_names_tex(2) = {'C'};
M_.endo_names_long(2) = {'C'};
M_.endo_names(3) = {'L'};
M_.endo_names_tex(3) = {'L'};
M_.endo_names_long(3) = {'L'};
M_.endo_names(4) = {'W'};
M_.endo_names_tex(4) = {'W'};
M_.endo_names_long(4) = {'W'};
M_.endo_names(5) = {'RK'};
M_.endo_names_tex(5) = {'RK'};
M_.endo_names_long(5) = {'RK'};
M_.endo_names(6) = {'K'};
M_.endo_names_tex(6) = {'K'};
M_.endo_names_long(6) = {'K'};
M_.endo_names(7) = {'I'};
M_.endo_names_tex(7) = {'I'};
M_.endo_names_long(7) = {'I'};
M_.endo_names(8) = {'LAMBDA'};
M_.endo_names_tex(8) = {'LAMBDA'};
M_.endo_names_long(8) = {'LAMBDA'};
M_.endo_names(9) = {'Q'};
M_.endo_names_tex(9) = {'Q'};
M_.endo_names_long(9) = {'Q'};
M_.endo_names(10) = {'MC'};
M_.endo_names_tex(10) = {'MC'};
M_.endo_names_long(10) = {'MC'};
M_.endo_names(11) = {'MT'};
M_.endo_names_tex(11) = {'MT'};
M_.endo_names_long(11) = {'MT'};
M_.endo_names(12) = {'PI'};
M_.endo_names_tex(12) = {'PI'};
M_.endo_names_long(12) = {'PI'};
M_.endo_names(13) = {'i'};
M_.endo_names_tex(13) = {'i'};
M_.endo_names_long(13) = {'i'};
M_.endo_names(14) = {'R'};
M_.endo_names_tex(14) = {'R'};
M_.endo_names_long(14) = {'R'};
M_.endo_names(15) = {'LRR'};
M_.endo_names_tex(15) = {'LRR'};
M_.endo_names_long(15) = {'LRR'};
M_.endo_partitions = struct();
M_.param_names = cell(10,1);
M_.param_names_tex = cell(10,1);
M_.param_names_long = cell(10,1);
M_.param_names(1) = {'BETA'};
M_.param_names_tex(1) = {'BETA'};
M_.param_names_long(1) = {'BETA'};
M_.param_names(2) = {'ETA'};
M_.param_names_tex(2) = {'ETA'};
M_.param_names_long(2) = {'ETA'};
M_.param_names(3) = {'DELTA'};
M_.param_names_tex(3) = {'DELTA'};
M_.param_names_long(3) = {'DELTA'};
M_.param_names(4) = {'ALPHA'};
M_.param_names_tex(4) = {'ALPHA'};
M_.param_names_long(4) = {'ALPHA'};
M_.param_names(5) = {'NU'};
M_.param_names_tex(5) = {'NU'};
M_.param_names_long(5) = {'NU'};
M_.param_names(6) = {'PSI'};
M_.param_names_tex(6) = {'PSI'};
M_.param_names_long(6) = {'PSI'};
M_.param_names(7) = {'RHOM'};
M_.param_names_tex(7) = {'RHOM'};
M_.param_names_long(7) = {'RHOM'};
M_.param_names(8) = {'HABIT'};
M_.param_names_tex(8) = {'HABIT'};
M_.param_names_long(8) = {'HABIT'};
M_.param_names(9) = {'KAPPA'};
M_.param_names_tex(9) = {'KAPPA'};
M_.param_names_long(9) = {'KAPPA'};
M_.param_names(10) = {'OMEGA'};
M_.param_names_tex(10) = {'OMEGA'};
M_.param_names_long(10) = {'OMEGA'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 15;
M_.param_nbr = 10;
M_.orig_endo_nbr = 15;
M_.aux_vars = [];
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
M_.orig_eq_nbr = 15;
M_.eq_nbr = 15;
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
 0 5 0;
 1 6 20;
 0 7 0;
 0 8 0;
 0 9 21;
 2 10 22;
 3 11 23;
 0 12 24;
 0 13 25;
 0 14 0;
 4 15 0;
 0 16 26;
 0 17 0;
 0 18 0;
 0 19 0;]';
M_.nstatic = 7;
M_.nfwrd   = 4;
M_.npred   = 1;
M_.nboth   = 3;
M_.nsfwrd   = 7;
M_.nspred   = 4;
M_.ndynamic   = 8;
M_.dynamic_tmp_nbr = [1; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'LAMBDA' ;
  4 , 'name' , 'L' ;
  5 , 'name' , 'RK' ;
  6 , 'name' , 'MC' ;
  7 , 'name' , 'PI' ;
  8 , 'name' , 'i' ;
  9 , 'name' , 'MT' ;
  10 , 'name' , 'Y' ;
  11 , 'name' , 'K' ;
  12 , 'name' , 'Q' ;
  13 , 'name' , '13' ;
  14 , 'name' , 'R' ;
  15 , 'name' , '15' ;
};
M_.mapping.Y.eqidx = [2 4 6 10 ];
M_.mapping.C.eqidx = [1 10 15 ];
M_.mapping.L.eqidx = [4 5 ];
M_.mapping.W.eqidx = [2 5 6 ];
M_.mapping.RK.eqidx = [5 13 ];
M_.mapping.K.eqidx = [2 4 5 6 11 13 ];
M_.mapping.I.eqidx = [10 11 12 ];
M_.mapping.LAMBDA.eqidx = [1 2 3 13 ];
M_.mapping.Q.eqidx = [12 13 ];
M_.mapping.MC.eqidx = [6 7 ];
M_.mapping.MT.eqidx = [8 9 ];
M_.mapping.PI.eqidx = [3 7 8 14 ];
M_.mapping.i.eqidx = [3 8 14 ];
M_.mapping.R.eqidx = [14 ];
M_.mapping.LRR.eqidx = [15 ];
M_.mapping.eps_m.eqidx = [9 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [2 6 7 11 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(15, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(10, 1);
M_.endo_trends = struct('deflator', cell(15, 1), 'log_deflator', cell(15, 1), 'growth_factor', cell(15, 1), 'log_growth_factor', cell(15, 1));
M_.NNZDerivatives = [55; 0; -1; ];
M_.static_tmp_nbr = [1; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load PARAMFILE;
set_param_value('BETA',BETA);
set_param_value('ETA',ETA);
set_param_value('NU',NU);
set_param_value('DELTA',DELTA);
set_param_value('ALPHA',ALPHA);
set_param_value('RHOM',RHOM);
set_param_value('PSI',PSI);
set_param_value('OMEGA',OMEGA);
set_param_value('HABIT',HABIT);
set_param_value('KAPPA',KAPPA);
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.01)^2;
options_.irf = 40;
options_.nograph = true;
options_.nomoments = true;
options_.order = 1;
var_list_ = {'K';'Y';'C';'PI';'MC';'i';'R';'I';'MT';'LRR';'RK'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('Brault_Khan_JMCB_2019_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Brault_Khan_JMCB_2019_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Brault_Khan_JMCB_2019_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Brault_Khan_JMCB_2019_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Brault_Khan_JMCB_2019_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('Brault_Khan_JMCB_2019_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('Brault_Khan_JMCB_2019_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
