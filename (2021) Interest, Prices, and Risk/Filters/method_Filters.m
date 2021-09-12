%% 1) Filters: HP, Baxter-King, Christiano-Fitzgerald

% Ex-ante real interest rate            
r_exAnte    = nrList - piiLRList;

% Apply filters
nSpecs = size(nrList,2);
datasetFilters = struct();
for iSpec = 1:nSpecs

    % Create timetable
    ts_Filters = timetable(ts.Properties.RowTimes, nrList(:,iSpec), piiLRList(:,iSpec), r_exAnte(:,iSpec), 'VariableNames', {'nr', 'pii', 'r'});
    ts_Filters = rmmissing(ts_Filters);

    % HP filter
    smoothingParam = 1600; % Yearly ? 100; Quarterly ? 1600; Monthly ? 14400
    [ts_Filters.trend_HP, ts_Filters.cycle_HP] = hpfilter(ts_Filters.r, smoothingParam);

    % Baxter-King filter
    % Use code from Pawe? Kowal (2005)
    % For the parameters: http://help.prognoz.com/en/mergedProjects/Lib/02_time_series_analysis/uimodelling_baxterkingfilter.htm
    %addpath(genpath(edu_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/filters', 'C:\')));
    BK_freqMin = 6;
    BK_freqMax = 32;
    BK_order   = 12;
    if length(ts_Filters.r) - BK_order >= BK_order + 1
        ts_Filters.cycle_BK = [NaN(2*BK_order,1); BK(ts_Filters.r, BK_freqMin, BK_freqMax, BK_order)];
        ts_Filters.trend_BK = ts_Filters.r - ts_Filters.cycle_BK;
    else
        ts_Filters.cycle_BK = NaN(length(ts_Filters.r),1);
        ts_Filters.trend_BK = NaN(length(ts_Filters.r),1);
    end

    % Christiano-Fitzgerald filter
    % Use code from Pawe? Kowal (2005)
    % For the parameters: http://help.prognoz.com/en/mergedProjects/Lib/02_time_series_analysis/uimodelling_baxterkingfilter.htm
    %addpath(genpath(edu_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/filters', 'C:\')));
    CF_freqMin = 6;
    CF_freqMax = 32;
    ts_Filters.cycle_CF = CF(ts_Filters.r, CF_freqMin, CF_freqMax);
    ts_Filters.trend_CF = ts_Filters.r - ts_Filters.cycle_CF;

    % Median estimation
    ts_Filters.cycle_Median = median([ts_Filters.cycle_HP, ts_Filters.cycle_BK, ts_Filters.cycle_CF], 2, 'omitnan');
    ts_Filters.trend_Median = median([ts_Filters.trend_HP, ts_Filters.trend_BK, ts_Filters.trend_CF], 2, 'omitnan');

    ts_Filters.cycle_Min = min([ts_Filters.cycle_HP, ts_Filters.cycle_BK, ts_Filters.cycle_CF], [], 2, 'omitnan');
    ts_Filters.trend_Min = min([ts_Filters.trend_HP, ts_Filters.trend_BK, ts_Filters.trend_CF], [], 2, 'omitnan');

    ts_Filters.cycle_Max = max([ts_Filters.cycle_HP, ts_Filters.cycle_BK, ts_Filters.cycle_CF], [], 2, 'omitnan');
    ts_Filters.trend_Max = max([ts_Filters.trend_HP, ts_Filters.trend_BK, ts_Filters.trend_CF], [], 2, 'omitnan');

    datasetFilters.(['Spec_' num2str(iSpec)]) = ts_Filters;
end