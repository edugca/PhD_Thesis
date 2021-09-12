function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = Brault_Khan_JMCB_2019.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(15, 27);
g1(1,1)=(-params(8));
g1(1,6)=1+params(1)*params(8);
g1(1,20)=(-(params(1)*params(8)));
g1(1,12)=(1-params(8))*(1-params(1)*params(8));
g1(2,5)=(-(params(2)/(1-params(4))));
g1(2,8)=1;
g1(2,2)=params(2)*params(4)/(1-params(4));
g1(2,12)=1;
g1(3,12)=1;
g1(3,24)=(-1);
g1(3,26)=1;
g1(3,17)=(-1);
g1(4,5)=(-(1/(1-params(4))));
g1(4,7)=1;
g1(4,2)=T(1);
g1(5,7)=(-(params(3)+1/params(1)-1));
g1(5,8)=(-(params(3)+1/params(1)-1));
g1(5,9)=1;
g1(5,2)=params(3)+1/params(1)-1;
g1(6,5)=(-T(1));
g1(6,8)=(-1);
g1(6,2)=T(1);
g1(6,14)=1;
g1(7,14)=(-params(6));
g1(7,16)=1;
g1(7,26)=(-params(1));
g1(8,15)=(-1);
g1(8,16)=(-params(5));
g1(8,17)=1;
g1(9,4)=(-params(7));
g1(9,15)=1;
g1(9,27)=(-1);
g1(10,5)=1;
g1(10,6)=(-(1-15*params(3)));
g1(10,11)=(-(15*params(3)));
g1(11,2)=(-(1-params(3)));
g1(11,10)=1;
g1(11,11)=(-params(3));
g1(12,3)=params(10);
g1(12,11)=(-(params(10)*(1+params(1))));
g1(12,23)=params(1)*params(10);
g1(12,13)=1;
g1(13,21)=(-1);
g1(13,2)=(-(15*params(9)));
g1(13,10)=15*params(9)+15*params(9);
g1(13,22)=(-(15*params(9)));
g1(13,12)=1;
g1(13,24)=(-1);
g1(13,13)=1;
g1(13,25)=(-(1-params(3)));
g1(14,26)=1;
g1(14,17)=(-1);
g1(14,18)=1;
g1(15,6)=(-1);
g1(15,19)=(-1);

end
