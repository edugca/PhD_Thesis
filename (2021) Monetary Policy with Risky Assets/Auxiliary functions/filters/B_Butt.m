function    B  = B_Butt(freq_min,freq_max,T,options)
% PURPOSE:  computes optimal filtering weights, B, in case of Butterworth
% filter using fourier transformation. Given time series X assuming
% infinite amount of data, filtered time series take the form
%
%       Y_t = sum _{j=-\infty}^\infty B_j X_{t+j}
%
% INPUT: in:
%        freq_min   - cut-off frequency - number of periods per cycle
%        freq_max   - cut-off frequency - number of periods per cycle
%        T          - number of computed weights.
%        options    - additional options
%                       n:      order of the Butterworth filter;
%                       int     options of integration algorithm
%                           type
%                               'D'     discretization
%                               'S'     adaptive Simpson quadrature
%                               'L'     adaptive Lobatto quadrature
%                               on default int = 'D'
%                           K   in case of discretization algorithm number
%                               of points, on default K = 10^4
%                           tol in case of Simpson and Lobatto quadrature
%                               absolute error tolerance, on default tol =
%                               10^-6
%
% OUTPUT: 
%        B      - optimal filtering weights, B=[B_0;l B_1; ... B_N-1]
%                 additionally B_{-k} = B_k
%
% ALGORITHM
%       optimal weight, B_j, is given by
%
%               B_j = 1/2*pi int_{-pi}^pi H(exp(-iw))*exp(iwj)dw
%
%       where H is a transfer function of the Butterworth filter.
%
% written by:
% Pawel Kowal
% Department of Economics
% Warsaw School of Economics
% pkowal3@sgh.waw.pl
       
n               = options.n;

B               = zeros(T,1);

freq_min        = max(freq_min,2); 
om_min          = 2*pi/freq_min;
om_max          = 2*pi/freq_max;
lambda_min      = tan(om_min/2)^(-2*n);
lambda_max      = tan(om_max/2)^(-2*n);
  
switch options.int.type
    case 'D'
        K           = options.int.K;

        s           = [-pi:2*pi/K:pi];
        s           = (s(1:end-1)+s(2:end))/2;
        p           = psi(s,n,lambda_min,lambda_max);

        for k=0:1:T-1
            B(k+1)  =real(sum(p.*exp(i*k*s))/2/pi*(2*pi/K));
        end   
        
    case 'S'
        tol         = options.int.tol;
        for k=0:1:T-1
            B(k+1)  =real(quad(@(x) psi(x,n,lambda_min,lambda_max,k),-pi,pi,tol)/2/pi);
        end  
    case 'L'
        tol         = options.int.tol;
        for k=0:1:T-1
            B(k+1)  =real(quadl(@(x) psi(x,n,lambda_min,lambda_max,k),-pi,pi,tol)/2/pi);
        end  
    otherwise
        error('unknown integration algorithm');
end

function out = psi(s,n,lambda_min,lambda_max,k)
s0              = exp(-i*s);
psi_min         = lambda_min.*(1-s0).^n.*(1-s0.^-1).^n;
psi_min         = psi_min./((1+s0).^n.*(1+s0.^-1).^n+psi_min);
psi_max         = lambda_max.*(1-s0).^n.*(1-s0.^-1).^n;
psi_max         = psi_max./((1+s0).^n.*(1+s0.^-1).^n+psi_max);
out             = psi_max-psi_min;

if nargin==5
    out         = out.*exp(i*k*s);
end