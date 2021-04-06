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
    T = smoothingCapital.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(10, 1);
lhs = (-y(6));
rhs = params(7)*y(2)-y(16)+(1-params(7))*params(6)*y(11)-y(17)+params(9)*y(13);
residual(1) = lhs - rhs;
lhs = (-y(6));
rhs = y(18)-y(16)+T(1)*(y(16)+(1+params(3))/(1-params(11))*y(15)-(1+params(11)*params(3))/(1-params(11))*y(7));
residual(2) = lhs - rhs;
lhs = y(11);
rhs = params(2)*y(17)+(y(6)+(params(11)+params(3))/(1-params(11))*y(5)-params(11)*(1+params(3))/(1-params(11))*y(1))*T(2);
residual(3) = lhs - rhs;
lhs = y(5);
rhs = y(7)*params(14)/params(13)+y(6)*(params(13)-params(14)*params(10))/params(13)-params(14)*(1-params(10))/params(13)*y(1);
residual(4) = lhs - rhs;
lhs = y(12);
rhs = params(12)*params(14)*(y(7)-y(1))-params(12)*params(14)*(y(1)-y(4));
residual(5) = lhs - rhs;
lhs = y(9);
rhs = params(9)*y(13)+params(7)*y(2)+(1-params(7))*params(6)*y(11);
residual(6) = lhs - rhs;
lhs = y(9);
rhs = y(17)+y(10);
residual(7) = lhs - rhs;
lhs = y(8);
rhs = y(4);
residual(8) = lhs - rhs;
lhs = y(13);
rhs = params(8)*y(3)+params(9)*x(it_, 1);
residual(9) = lhs - rhs;
lhs = y(14);
rhs = y(1);
residual(10) = lhs - rhs;

end
