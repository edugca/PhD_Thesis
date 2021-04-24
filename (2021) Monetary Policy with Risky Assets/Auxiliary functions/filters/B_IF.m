function  B  = B_IF(freq_min,freq_max,N)
% PURPOSE:  computes optimal filtering weights, B, in case of the ideal
% bandpass filter using fourier transformation. Given time series X assuming
% infinite amount of data, filtered time series take the form
%
%       Y_t = sum _{j=-\infty}^\infty B_j X_{t+j}
%
% INPUT: in:
%        freq_min   - cut-off frequency - minimal number of periods per cycle
%        freq_max   - cut-off frequency - maximal number of periods per
%                       cycle
%        N          - number of computed weights.
%
% OUTPUT: 
%        B      - optimal filtering weights, B=[B_0;l B_1; ... B_N-1]
%                 additionally B_{-k} = B_k
%
% ALGORITHM
%       optimal weight, B_j, is given by
%
%               B_0     = (b-a)/pi
%               B_j     = (sin(bj)-sin(aj)/(pi x j)
%
%       where [a,b] is an interval of frequencies that are passed by the
%       filter
%
% written by:
% Pawel Kowal
% Department of Economics
% Warsaw School of Economics
% pkowal3@sgh.waw.pl
            
freq_min    = max(2,freq_min);
a           = 2*pi/freq_max;
b           = 2*pi/freq_min;

J           = [1:1:N-1]';
B           = (sin(b*J)-sin(a*J))./(pi*J);
B           = [(b-a)/pi;B];

