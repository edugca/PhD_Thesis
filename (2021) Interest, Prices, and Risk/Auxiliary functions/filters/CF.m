function  Y  = CF(X,freq_min,freq_max)
% PURPOSE:  computes cyclical component of given time series using the
% Christiano-Fitzgerald bandpass filter
%
% INPUT: 
%        X          - time series (T x K)
%        freq_min   - cut-off frequency - number of periods per cycle
%        freq_max   - cut-off frequency - number of periods per cycle
%
% OUTPUT: 
%        Y      - cyclical component of X
%
% written by:
% Pawel Kowal
% Department of Economics
% Warsaw School of Economics
% pkowal3@sgh.waw.pl
           
freq_min        = max(2,freq_min);

a               = 2*pi/freq_max;
b               = 2*pi/freq_min;

[T,nx]          = size(X);
% remove ols trend
dX              = (X(end,:)-X(1,:))/(T-1);
X               = X - [0:1:T-1]'*dX;

J               = [1:1:T]';
B               = (sin(b*J)-sin(a*J))./(pi*J);
B               = [(b-a)/pi;B];

B_main          = B(1:T);
B_tmp           = B(1:T);

n2              = T-1;
n1              = 0;
Y               = 0*X;
for i=1:1:T   
    B_end       = sum(B_tmp(1+n2:end));
    B_start     = sum(B_tmp(1+n1:end));

    BB          = B_main;
    BB(1)       = B_start;
    BB(end)     = B_end;

    Y(i,:)      = sum(X.*(BB*ones(1,nx)));

    B_main(end) = [];
    B_main      =[B(i+1);B_main];

    n2          = n2-1;
    n1          = n1+1;    
end