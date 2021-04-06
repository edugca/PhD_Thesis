function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = risky_ClosedEconomy.static_resid_tt(T, y, x, params);
end
residual = zeros(16, 1);
lhs = y(4);
rhs = params(2)*y(4)+T(1)*y(1);
residual(1) = lhs - rhs;
lhs = y(1);
rhs = y(1)-((1-params(2)*params(13)*params(14))*y(2)-y(4)+(params(2)*params(13)*params(14)-params(13))/(1-params(13))*y(9)+params(2)*params(13)*params(14)*y(10));
residual(2) = lhs - rhs;
lhs = y(1);
rhs = y(1)-(y(6)-y(4));
residual(3) = lhs - rhs;
lhs = y(2);
rhs = y(11)+y(2)*params(20)+(1-params(20))*(params(6)*y(4)+y(1)*params(19));
residual(4) = lhs - rhs;
residual(5) = y(5);
lhs = y(3);
rhs = (1-params(13))*params(15)/(1+params(15))*y(2)-y(4)+params(13)*(params(14)+(-1)-params(15))/(1+params(15))*y(9)+params(13)*params(14)/(1+params(15))*y(10);
residual(6) = lhs - rhs;
lhs = y(7);
rhs = y(6)-y(4);
residual(7) = lhs - rhs;
lhs = y(8);
rhs = y(2)-y(6);
residual(8) = lhs - rhs;
lhs = y(9);
rhs = y(12);
residual(9) = lhs - rhs;
lhs = y(10);
rhs = y(13);
residual(10) = lhs - rhs;
lhs = y(11);
rhs = params(7)*y(11)+params(10)*x(1);
residual(11) = lhs - rhs;
lhs = y(12);
rhs = params(8)*y(12)+params(11)*x(2);
residual(12) = lhs - rhs;
lhs = y(13);
rhs = params(9)*y(13)+params(12)*x(3);
residual(13) = lhs - rhs;
lhs = y(14);
rhs = y(11);
residual(14) = lhs - rhs;
lhs = y(15);
rhs = y(11)*params(6)*T(2);
residual(15) = lhs - rhs;
lhs = y(16);
rhs = y(11)*T(2)*(-params(7));
residual(16) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
