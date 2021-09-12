function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = Brault_Khan_JMCB_2019.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(15, 1);
lhs = (1-params(8))*(1-params(1)*params(8))*y(12);
rhs = params(8)*y(1)-(1+params(1)*params(8))*y(6)+params(1)*params(8)*y(20);
residual(1) = lhs - rhs;
lhs = y(12)+y(8);
rhs = params(2)/(1-params(4))*y(5)-params(2)*params(4)/(1-params(4))*y(2);
residual(2) = lhs - rhs;
lhs = y(12);
rhs = y(24)+y(17)-y(26);
residual(3) = lhs - rhs;
lhs = y(7);
rhs = y(5)*1/(1-params(4))-y(2)*T(1);
residual(4) = lhs - rhs;
lhs = y(9);
rhs = (params(3)+1/params(1)-1)*(y(7)+y(8)-y(2));
residual(5) = lhs - rhs;
lhs = y(14);
rhs = y(8)+y(5)*T(1)-y(2)*T(1);
residual(6) = lhs - rhs;
lhs = y(16);
rhs = params(1)*y(26)+y(14)*params(6);
residual(7) = lhs - rhs;
lhs = y(17);
rhs = y(16)*params(5)+y(15);
residual(8) = lhs - rhs;
lhs = y(15);
rhs = params(7)*y(4)+x(it_, 1);
residual(9) = lhs - rhs;
lhs = y(5);
rhs = y(6)*(1-15*params(3))+y(11)*15*params(3);
residual(10) = lhs - rhs;
lhs = y(10);
rhs = y(2)*(1-params(3))+params(3)*y(11);
residual(11) = lhs - rhs;
lhs = y(13);
rhs = y(11)*params(10)*(1+params(1))-params(10)*y(3)-params(1)*params(10)*y(23);
residual(12) = lhs - rhs;
lhs = y(13)+(y(10)-y(2))*15*params(9);
rhs = y(24)-y(12)+y(21)+(1-params(3))*y(25)+(y(22)-y(10))*15*params(9);
residual(13) = lhs - rhs;
lhs = y(18);
rhs = y(17)-y(26);
residual(14) = lhs - rhs;
lhs = (-y(6));
rhs = y(19);
residual(15) = lhs - rhs;

end
