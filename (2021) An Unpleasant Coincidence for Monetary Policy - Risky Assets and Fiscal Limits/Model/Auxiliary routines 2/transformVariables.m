function [obsData, filtData, smoothData] = transformVariables(varName, mapEstimParams, obsData, filtData, smoothData)

if ismember(varName, {'c', 'y', 'g'})
    obsData     = obsData * 100;
    filtData    = (log(filtData./lagmatrix(filtData,1)) + edu_IIF(mapEstimParams.isKey('ytrend'), @() mapEstimParams('ytrend'), 0)) * 100;
    smoothData  = (log(smoothData./lagmatrix(smoothData,1)) + edu_IIF(mapEstimParams.isKey('ytrend'), @() mapEstimParams('ytrend'), 0)) * 100 ;
elseif ismember(varName, {'n'})
    obsData     = obsData * 100;
    filtData    = filtData/mapEstimParams('unempSS') * mean(obsData, 'omitnan');
    smoothData  = smoothData/mapEstimParams('unempSS') * mean(obsData, 'omitnan');
elseif ismember(varName, {'Pii'})
    obsData     = obsData * 100;
    filtData    = log(filtData) * 100;
    smoothData  = log(smoothData) * 100;
elseif ismember(varName, {'nrPolicy'})
    obsData     = obsData * 100;
    filtData    = ((1+filtData).^4 - 1) .* 100;
    smoothData  = ((1+smoothData).^4 - 1) .* 100;
elseif ismember(varName, {'rPolicy'})
    obsData     =  obsData .* 100;
    filtData    = ((1+filtData).^4 - 1) .* 100;
    smoothData  = ((1+smoothData).^4 - 1) .* 100;
elseif ismember(varName, {'tau', 'bY'})
    obsData     =  obsData ;
    filtData    = filtData .* 100;
    smoothData  = smoothData .* 100;
elseif ismember(varName, {'polDef'})
    obsData     = obsData * 100;
    filtData    = ( (filtData - lagmatrix(filtData,1)) + edu_IIF(mapEstimParams.isKey('polDefWedge'), @() mapEstimParams('polDefWedge'), 0) ) * 100;
    smoothData  = ( (smoothData - lagmatrix(smoothData,1)) + edu_IIF(mapEstimParams.isKey('polDefWedge'), @() mapEstimParams('polDefWedge'), 0) ) * 100;
end

end