function Y = butterworth(X,freq,n)
% PURPOSE:  computes the cyclical component of given time series using the Butterworth
% filter
%
% INPUT: 
%        X      - time series (T x K)
%        n      - order of the Butterworth filter
%        freq   - cut-off frequency - number of periods per cycle
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

cut_off         = 2*pi/freq;
mu              = (1/tan(cut_off/2))^(2*n);

I               = speye(T);
LT              = spdiags(ones(T-1,1),-1,T,T);
LT              = (I-LT)^n;
Q               = LT(3:end,:)';
SIGMA_R         = Q'*Q;
SIGMA_T         = abs(SIGMA_R);

g               = Q'*X;
b               = (SIGMA_T+mu*SIGMA_R)\g;
Y               = mu*Q*b;