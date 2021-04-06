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
    T = smoothingCapital.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(10, 19);
g1(1,6)=(-1);
g1(1,16)=1;
g1(1,2)=(-params(7));
g1(1,11)=(-((1-params(7))*params(6)));
g1(1,17)=1;
g1(1,13)=(-params(9));
g1(2,15)=(-(T(1)*(1+params(3))/(1-params(11))));
g1(2,6)=(-1);
g1(2,16)=(-(T(1)-1));
g1(2,7)=(-(T(1)*(-((1+params(11)*params(3))/(1-params(11))))));
g1(2,18)=(-1);
g1(3,5)=(-((params(11)+params(3))/(1-params(11))*T(2)));
g1(3,6)=(-T(2));
g1(3,1)=(-(T(2)*(-(params(11)*(1+params(3))/(1-params(11))))));
g1(3,11)=1;
g1(3,17)=(-params(2));
g1(4,5)=1;
g1(4,6)=(-((params(13)-params(14)*params(10))/params(13)));
g1(4,1)=params(14)*(1-params(10))/params(13);
g1(4,7)=(-(params(14)/params(13)));
g1(5,1)=(-((-(params(12)*params(14)))-params(12)*params(14)));
g1(5,7)=(-(params(12)*params(14)));
g1(5,12)=1;
g1(5,4)=(-(params(12)*params(14)));
g1(6,2)=(-params(7));
g1(6,9)=1;
g1(6,11)=(-((1-params(7))*params(6)));
g1(6,13)=(-params(9));
g1(7,9)=1;
g1(7,10)=(-1);
g1(7,17)=(-1);
g1(8,8)=1;
g1(8,4)=(-1);
g1(9,3)=(-params(8));
g1(9,13)=1;
g1(9,19)=(-params(9));
g1(10,1)=(-1);
g1(10,14)=1;

end
