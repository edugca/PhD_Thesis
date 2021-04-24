function yForecast = pub_ModelRollingForecast(estMdl, nPeriodsAhead, y)

nObs = length(y);
yForecast = NaN(nObs, nPeriodsAhead);

% First observation
startObs = find(~isnan(y), 1) + estMdl.P;

for iObs = startObs:nObs
    yForecast(iObs,:) = forecast(estMdl,nPeriodsAhead,y(1:iObs-1));
end

end