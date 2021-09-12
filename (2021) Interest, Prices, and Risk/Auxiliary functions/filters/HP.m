function Y = HP(X,freq,type)
% PURPOSE:  computes cyclical component of given time series using the HP
% filter
%
% INPUT: 
%        X      - time series (T x n)
%        type   - filter type
%                   0   - cyclical component is determined by the parameter
%                         lambda (default)
%                   1   - cyclical component is determined by the cut-off
%                         frequency
%        freq   - if type = 0 then freq is the lambda parameter of the HP filter
%                           suggested values
%                           6       for annual data
%                           1600    for quarterly data
%                           129600  for monthly data
%                 if type = 1 then freq is the cut-off frequency (number of
%                       periods per cycle) of the filter. Cut-off frequency is 
%                       a frequency when filter gain is equal 0.5.
%                       Dependence between the parameter lambda and the cut-off
%                       frequency:
%                               lambda = (2*sin(pi/freq))^-4
%                               freq   = pi/arcsin(1/2*lambda^-1/4)
%
% OUTPUT: 
%        Y      - cyclical component of X
%
% written by:
% Pawel Kowal
% Department of Economics
% Warsaw School of Economics
% pkowal3@sgh.waw.pl

T               = size(X,1);

I               = speye(T);
LT              = spdiags(ones(T-1,1),-1,T,T);
LT              = (I-LT)^2;
Q               = LT(3:end,:)';
SIGMA_R         = Q'*Q;
SIGMA_T         = speye(T-2);

lambda          = freq;
if nargin>2 && type ==1
    lambda      = (2*sin(pi/freq))^-4;
end
    
g               = Q'*X;
b               = (SIGMA_T+lambda*SIGMA_R)\g;
Y               = lambda*Q*b;