function  Y  = BK(X,freq_min,freq_max,K)
% PURPOSE:  computes cyclical component of given time series using the
% Baxter-King bandpass filter
%
% INPUT: 
%        X          - time series (T x K)
%        freq_min   - cut-off frequency - number of periods per cycle
%        freq_max   - cut-off frequency - number of periods per cycle
%        K          - order of the filter
%
% OUTPUT: 
%        Y      - cyclical component of X
%
% written by:
% Pawel Kowal
% Department of Economics
% Warsaw School of Economics
% pkowal3@sgh.waw.pl
           
freq_min    = max(2,freq_min);

a           = 2*pi/freq_max;
b           = 2*pi/freq_min;

J           = [1:1:K]';
B           = (sin(b*J)-sin(a*J))./(pi*J);
B           = [(b-a)/pi;B];

%Baxter-King approximation
s           = B(1) + 2*sum(B(2:end));
BB          = B - s/(2*K+1);
BB          = [flipud(BB(2:end));BB];

%filtering
T           = size(X,1);
for t = K+1:1:T-K
    Y(t,:)  = BB'*X(t-K:t+K,:);    
end
Y(1:K,:)=[];