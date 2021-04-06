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
    T = risky_ClosedEconomy.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(16, 25);
g1(1,5)=(-T(1));
g1(1,8)=1;
g1(1,22)=(-params(2));
g1(2,5)=1;
g1(2,21)=(-1);
g1(2,6)=1-params(2)*params(13)*params(14);
g1(2,22)=(-1);
g1(2,13)=(params(2)*params(13)*params(14)-params(13))/(1-params(13));
g1(2,14)=params(2)*params(13)*params(14);
g1(3,5)=1;
g1(3,21)=(-1);
g1(3,22)=(-1);
g1(3,10)=1;
g1(4,5)=(-((1-params(20))*params(19)));
g1(4,1)=(-params(20));
g1(4,6)=1;
g1(4,8)=(-(params(6)*(1-params(20))));
g1(4,15)=(-1);
g1(5,9)=1;
g1(6,6)=(-((1-params(13))*params(15)/(1+params(15))));
g1(6,7)=1;
g1(6,22)=1;
g1(6,13)=(-(params(13)*(params(14)+(-1)-params(15))/(1+params(15))));
g1(6,14)=(-(params(13)*params(14)/(1+params(15))));
g1(7,22)=1;
g1(7,10)=(-1);
g1(7,11)=1;
g1(8,6)=(-1);
g1(8,10)=1;
g1(8,12)=1;
g1(9,13)=1;
g1(9,16)=(-1);
g1(10,14)=1;
g1(10,17)=(-1);
g1(11,2)=(-params(7));
g1(11,15)=1;
g1(11,23)=(-params(10));
g1(12,3)=(-params(8));
g1(12,16)=1;
g1(12,24)=(-params(11));
g1(13,4)=(-params(9));
g1(13,17)=1;
g1(13,25)=(-params(12));
g1(14,15)=(-1);
g1(14,18)=1;
g1(15,15)=(-(params(6)*T(2)));
g1(15,19)=1;
g1(16,15)=(-((-params(7))*T(2)));
g1(16,20)=1;

end
