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
    T = Brault_Khan_JMCB_2019.static_g1_tt(T, y, x, params);
end
g1 = zeros(15, 15);
g1(1,2)=(-(params(1)*params(8)+params(8)-(1+params(1)*params(8))));
g1(1,8)=(1-params(8))*(1-params(1)*params(8));
g1(2,1)=(-(params(2)/(1-params(4))));
g1(2,4)=1;
g1(2,6)=params(2)*params(4)/(1-params(4));
g1(2,8)=1;
g1(3,12)=1;
g1(3,13)=(-1);
g1(4,1)=(-(1/(1-params(4))));
g1(4,3)=1;
g1(4,6)=T(1);
g1(5,3)=(-(params(3)+1/params(1)-1));
g1(5,4)=(-(params(3)+1/params(1)-1));
g1(5,5)=1;
g1(5,6)=params(3)+1/params(1)-1;
g1(6,1)=(-T(1));
g1(6,4)=(-1);
g1(6,6)=T(1);
g1(6,10)=1;
g1(7,10)=(-params(6));
g1(7,12)=1-params(1);
g1(8,11)=(-1);
g1(8,12)=(-params(5));
g1(8,13)=1;
g1(9,11)=1-params(7);
g1(10,1)=1;
g1(10,2)=(-(1-15*params(3)));
g1(10,7)=(-(15*params(3)));
g1(11,6)=1-(1-params(3));
g1(11,7)=(-params(3));
g1(12,7)=(-(params(10)*(1+params(1))-params(10)-params(1)*params(10)));
g1(12,9)=1;
g1(13,5)=(-1);
g1(13,9)=1-(1-params(3));
g1(14,12)=1;
g1(14,13)=(-1);
g1(14,14)=1;
g1(15,2)=(-1);
g1(15,15)=(-1);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
