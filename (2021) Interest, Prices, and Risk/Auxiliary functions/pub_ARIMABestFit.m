function bestMdl = pub_ARIMABestFit(y, pMax, dMax, qMax, infoCriterion, ...
    seasonality, sarMax, smaMax)

% INPUT
% infoCriterion = 'aic' or 'bic'
%
% OUTPUT
% bestMdl       = [p, d, q]

LOGL    = zeros(pMax + 1, dMax + 1, qMax + 1); % Initialize
PQ      = zeros(pMax + 1, dMax + 1, qMax + 1); % Initialize
for p = 0:pMax
    for d = 0:dMax
       for q = 0:qMax
           
            mdl = arima(p, d, q);
            mdl.Seasonality = seasonality;
            %mdl.SAR = sarMax;
            %mdl.SMA = smaMax;
            
            try
                [~,~,logL]    = estimate(mdl,y,'Display','off');
            catch
               logL = nan; 
            end
            
            LOGL(p + 1, d + 1, q + 1)     = logL;
            PQ(p + 1, d + 1, q + 1)       = p+q;
            
       end
    end
end

LOGL    = reshape(LOGL,[],1);
PQ      = reshape(PQ,[],1);

nObs = sum(~isnan(y));
[aic,bic] = aicbic(LOGL, PQ+1 + seasonality, nObs);

if strcmp(infoCriterion, 'aic')
    aic = reshape(aic, pMax + 1, dMax + 1,qMax + 1);
    
    [~, I] = min(aic(:));
    [I1,I2,I3] = ind2sub(size(aic),I);
    
    bestMdl = [I1-1, I2-1, I3-1];
elseif strcmp(infoCriterion, 'bic')
    bic = reshape(bic, pMax + 1, dMax + 1,qMax + 1);
    
    [~, I] = min(bic(:));
    [I1,I2,I3] = ind2sub(size(bic),I);
    
    bestMdl = [I1-1, I2-1, I3-1];
end

end