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
    T = smoothingCapital.static_g1_tt(T, y, x, params);
end
g1 = zeros(10, 10);
g1(1,5)=(-params(7));
g1(1,7)=(-((1-params(7))*params(6)-1));
g1(1,9)=(-params(9));
g1(2,1)=(-(T(1)*(1+params(3))/(1-params(11))));
g1(2,2)=(-1)-(T(1)-1);
g1(2,3)=(-(T(1)*(-((1+params(11)*params(3))/(1-params(11))))));
g1(2,8)=(-1);
g1(3,1)=(-(T(2)*(params(11)+params(3))/(1-params(11))));
g1(3,2)=(-T(2));
g1(3,3)=(-(T(2)*(-(params(11)*(1+params(3))/(1-params(11))))));
g1(3,7)=1-params(2);
g1(4,1)=1;
g1(4,2)=(-((params(13)-params(14)*params(10))/params(13)));
g1(4,3)=(-(params(14)/params(13)-params(14)*(1-params(10))/params(13)));
g1(5,8)=1;
g1(6,5)=1-params(7);
g1(6,7)=(-((1-params(7))*params(6)));
g1(6,9)=(-params(9));
g1(7,5)=1;
g1(7,6)=(-1);
g1(7,7)=(-1);
g1(8,3)=(-1);
g1(8,4)=1;
g1(9,9)=1-params(8);
g1(10,3)=(-1);
g1(10,10)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
