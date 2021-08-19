function [p,priors] = paramsForEstimation_MLE(start_at_mode)
if nargin<1
    start_at_mode=false;
end

p = struct();

%p.beta = 0.99;
%p.mu = 0.34;
%p.gam = 2;

priors=struct();

if start_at_mode
    priors.sigma = {1,  0.1,    5};
else
    priors.sigma = {1,  0.1,    5};
end

% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);
for ip=1:numel(fields)
    name=fields{ip};
    p.(name)=priors.(name){1};
end

end