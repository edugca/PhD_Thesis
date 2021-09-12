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
    T = Brault_Khan_JMCB_2019.static_resid_tt(T, y, x, params);
end
residual = zeros(15, 1);
lhs = (1-params(8))*(1-params(1)*params(8))*y(8);
rhs = params(8)*y(2)-y(2)*(1+params(1)*params(8))+params(1)*params(8)*y(2);
residual(1) = lhs - rhs;
lhs = y(8)+y(4);
rhs = params(2)/(1-params(4))*y(1)-params(2)*params(4)/(1-params(4))*y(6);
residual(2) = lhs - rhs;
lhs = y(8);
rhs = y(8)+y(13)-y(12);
residual(3) = lhs - rhs;
lhs = y(3);
rhs = y(1)*1/(1-params(4))-y(6)*T(1);
residual(4) = lhs - rhs;
lhs = y(5);
rhs = (params(3)+1/params(1)-1)*(y(3)+y(4)-y(6));
residual(5) = lhs - rhs;
lhs = y(10);
rhs = y(4)+y(1)*T(1)-y(6)*T(1);
residual(6) = lhs - rhs;
lhs = y(12);
rhs = params(1)*y(12)+y(10)*params(6);
residual(7) = lhs - rhs;
lhs = y(13);
rhs = y(12)*params(5)+y(11);
residual(8) = lhs - rhs;
lhs = y(11);
rhs = y(11)*params(7)+x(1);
residual(9) = lhs - rhs;
lhs = y(1);
rhs = (1-15*params(3))*y(2)+15*params(3)*y(7);
residual(10) = lhs - rhs;
lhs = y(6);
rhs = y(6)*(1-params(3))+params(3)*y(7);
residual(11) = lhs - rhs;
lhs = y(9);
rhs = y(7)*params(10)*(1+params(1))-y(7)*params(10)-y(7)*params(1)*params(10);
residual(12) = lhs - rhs;
lhs = y(9);
rhs = y(5)+(1-params(3))*y(9);
residual(13) = lhs - rhs;
lhs = y(14);
rhs = y(13)-y(12);
residual(14) = lhs - rhs;
lhs = (-y(2));
rhs = y(15);
residual(15) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
