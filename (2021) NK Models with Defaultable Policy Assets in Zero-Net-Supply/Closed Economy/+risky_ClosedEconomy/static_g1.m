function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = risky_ClosedEconomy.static_g1_tt(T, y, x, params);
end
g1 = zeros(16, 16);
g1(1,1)=(-T(1));
g1(1,4)=1-params(2);
g1(2,2)=1-params(2)*params(13)*params(14);
g1(2,4)=(-1);
g1(2,9)=(params(2)*params(13)*params(14)-params(13))/(1-params(13));
g1(2,10)=params(2)*params(13)*params(14);
g1(3,4)=(-1);
g1(3,6)=1;
g1(4,1)=(-((1-params(20))*params(19)));
g1(4,2)=1-params(20);
g1(4,4)=(-(params(6)*(1-params(20))));
g1(4,11)=(-1);
g1(5,5)=1;
g1(6,2)=(-((1-params(13))*params(15)/(1+params(15))));
g1(6,3)=1;
g1(6,4)=1;
g1(6,9)=(-(params(13)*(params(14)+(-1)-params(15))/(1+params(15))));
g1(6,10)=(-(params(13)*params(14)/(1+params(15))));
g1(7,4)=1;
g1(7,6)=(-1);
g1(7,7)=1;
g1(8,2)=(-1);
g1(8,6)=1;
g1(8,8)=1;
g1(9,9)=1;
g1(9,12)=(-1);
g1(10,10)=1;
g1(10,13)=(-1);
g1(11,11)=1-params(7);
g1(12,12)=1-params(8);
g1(13,13)=1-params(9);
g1(14,11)=(-1);
g1(14,14)=1;
g1(15,11)=(-(params(6)*T(2)));
g1(15,15)=1;
g1(16,11)=(-(T(2)*(-params(7))));
g1(16,16)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
