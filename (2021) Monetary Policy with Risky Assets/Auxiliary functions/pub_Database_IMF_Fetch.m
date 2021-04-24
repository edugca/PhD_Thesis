function data = pub_Database_IMF_Fetch(database, indicator, country2L, frequency, startDate, endDate)

% database  = 'IFS'
% indicator = 'FPOLM_PA'
% country2L = 'US'
% frequency = 'A' is annual, 'Q' is quarterly, 'M' is monthly
% startDate = '2000'
% endDate   = '2010'

% Load IMF routine
pathRoutine = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Robert Kirkby/IMF/', 'C:\');
addpath(pathRoutine);

structIMF = getIMFData(...
    database, indicator, ...
    country2L, frequency, ...
    startDate, endDate);

% Avoid error when this code is executed sequentially
pause(1)

% Check whether it was successful
if sum(isnan(structIMF.Data)) ~= 2

    % Check if returns is a cell vector or a struct
    if iscell(structIMF.IMFcodes.Obs)
        dates = cellfun(@(x) x.x_TIME_PERIOD, structIMF.IMFcodes.Obs, 'UniformOutput', false);
    else
        dates = {structIMF.IMFcodes.Obs.x_TIME_PERIOD}';
    end
    
    % Convert dates format
    if strcmp(frequency, 'M')
        dates = datetime(dates, 'InputFormat', 'yyyy-MM');
    elseif strcmp(frequency, 'Q')
        dates = datetime(dates, 'InputFormat', 'yyyy-QQQ');
    end

    % The actual dates and data are kept in:
    data = timetable(dates, structIMF.Data(:,2), 'VariableNames', {indicator});
   
    disp(['Fetched from IMF: ' database ' ' indicator ' ' country2L]);

else
    
   data = nan; 
   
   disp(['Failed to fetch from IMF: ' database ' ' indicator ' ' country2L]);
   
end

end