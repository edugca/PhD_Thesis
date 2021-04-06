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
    T = risky_ClosedEconomy.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(16, 1);
lhs = y(8);
rhs = params(2)*y(22)+T(1)*y(5);
residual(1) = lhs - rhs;
lhs = y(5);
rhs = y(21)-((1-params(2)*params(13)*params(14))*y(6)-y(22)+(params(2)*params(13)*params(14)-params(13))/(1-params(13))*y(13)+params(2)*params(13)*params(14)*y(14));
residual(2) = lhs - rhs;
lhs = y(5);
rhs = y(21)-(y(10)-y(22));
residual(3) = lhs - rhs;
lhs = y(6);
rhs = y(15)+params(20)*y(1)+(1-params(20))*(params(6)*y(8)+y(5)*params(19));
residual(4) = lhs - rhs;
residual(5) = y(9);
lhs = y(7);
rhs = (1-params(13))*params(15)/(1+params(15))*y(6)-y(22)+params(13)*(params(14)+(-1)-params(15))/(1+params(15))*y(13)+params(13)*params(14)/(1+params(15))*y(14);
residual(6) = lhs - rhs;
lhs = y(11);
rhs = y(10)-y(22);
residual(7) = lhs - rhs;
lhs = y(12);
rhs = y(6)-y(10);
residual(8) = lhs - rhs;
lhs = y(13);
rhs = y(16);
residual(9) = lhs - rhs;
lhs = y(14);
rhs = y(17);
residual(10) = lhs - rhs;
lhs = y(15);
rhs = params(7)*y(2)+params(10)*x(it_, 1);
residual(11) = lhs - rhs;
lhs = y(16);
rhs = params(8)*y(3)+params(11)*x(it_, 2);
residual(12) = lhs - rhs;
lhs = y(17);
rhs = params(9)*y(4)+params(12)*x(it_, 3);
residual(13) = lhs - rhs;
lhs = y(18);
rhs = y(15);
residual(14) = lhs - rhs;
lhs = y(19);
rhs = y(15)*params(6)*T(2);
residual(15) = lhs - rhs;
lhs = y(20);
rhs = y(15)*(-params(7))*T(2);
residual(16) = lhs - rhs;

end
