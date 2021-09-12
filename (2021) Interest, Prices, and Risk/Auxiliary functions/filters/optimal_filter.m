function Y = optimal_filter(X,type,freq_min,freq_max,options)
% PURPOSE:  computes cyclical component of X using optimal finite sample
% approximation of the ideal filter.
% Given time series X assuming infinite amount of data, filtered time
% series take the form 
%
%       Y_t  = sum _{j=-\infty}^\infty B_j X_{t+j}
%
% where B_j depends on filter type. Finite sample approximation takes the
% form 
%
%       YY_t = sum _{j=-n1}^n2 BB_tj X_{t+j}
%
% where coefficients BB_tj are chosen to minimize mean square error for
% each t
%
%       min E{(Y_t-YY_t)^2}
%
% INPUT:
%        X          - time series (T x n)
%        type       - filter type
%                       'IF'    - ideal bandpass filter
%                       'B'     - Butterworth filter
%                       'HP'    - HP filter
%        freq_min   - cut-off frequency - number of periods per cycle
%        freq_max   - cut-off frequency - number of periods per cycle
%        options    - additional options
%                       type:   'S' for stationary and symmetric filter, 
%                               'A' for nonstationary and asymmetric filter (default)  
%                       nT:     number of nonzero coefficients BB_j, on default nT = T+1.
%                       nK:     for symmetric filters only, order of moving
%                               average representation of the filter,
%                               (2*K+1) is a number of nonzero coefficients
%                               BB_j, on default nK = 12.
%                       n:      for Butterworth filter only, filter order,
%                               on default n = 2.
%                       int     for Butterworth and HP filter, options of integration
%                               algorithm
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
%        Y          - cyclical component of X
%
% written by:
% Pawel Kowal
% Department of Economics
% Warsaw School of Economics
% pkowal3@sgh.waw.pl


%default options
default_options.type        = 'A';
default_options.nK          = 12;
default_options.n           = 2;
default_options.int.type    = 'D';
default_options.int.K       = 10^4+1;
default_options.int.tol     = 10^-6;


if nargin<5
    options         = default_options;
end

[T,nx]              = size(X);
B                   = zeros(T+1,1);

if isfield(options,'nT')
    nT              = options.nT;
else
    nT              = T+1;
end
if ~isfield(options,'n')
    options.n       = default_options.n;
end
if ~isfield(options,'int')
    options.int     = default_options.int;
end
if ~isfield(options.int,'type')
    options.int.type= default_options.int.type;
end
if ~isfield(options.int,'K')
    options.int.K   = default_options.int.K;
end
if ~isfield(options.int,'tol')
    options.int.tol = default_options.int.tol;
end

switch type
    case 'IF'
        B(1:nT)     = B_IF(freq_min,freq_max,nT); 
    case 'B'
        B(1:nT)     = B_Butt(freq_min,freq_max,nT,options); 
    case 'HP'
        B(1:nT)     = B_HP(freq_min,freq_max,nT,options); 
    otherwise
        error('unknown filter type');
end

Y                   = 0*X;
% remove ols trend
dX                  = (X(end,:)-X(1,:))/(T-1);
X                   = X - [0:1:T-1]'*dX;

if options.type=='A'
    %   ASYMMETRIC AND NONSTATIONARY FILTER
    B_main          = B(1:T);
    B_tmp           = B(1:T);

    n2              = T-1;
    n1              = 0;
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
else
    %   SYMMETRIC AND STATIONARY FILTER, BUT DATA LOSS
    if isfield(options,'nK')
        K           = options.nK;
    else
        K           = default_options.nK;
    end
    
    B               = B(1:K+1);
    s               = B(1) + 2*sum(B(2:end));
    B               = B - s/(2*K+1);
    B               = [flipud(B(2:end));B];

    for t = K+1:1:T-K
        Y(t,:)      = B'*X(t-K:t+K,:);    
    end
    Y(T-K+1:T,:)    =[];  
    Y(1:K,:)        =[];        
end