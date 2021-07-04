%% Fetch data from the SGS

function ts_final = pub_Database_SGS_Get(indicatorCode, startDate, endDate)

beginUrl = 'http://api.bcb.gov.br/dados/serie/bcdata.sgs.';

formatOut = 'dd/mm/yyyy';

if isempty(startDate)
    startDate_str = '01/01/1900';
else
    startDate_str = datestr(startDate,formatOut);
end

if isempty(endDate)
    endDate_str = '01/01/2100';
else
    endDate_str = datestr(endDate,formatOut);
end

if isnumeric(indicatorCode)
    indicatorCode = {indicatorCode};
end

ts_final = timetable();
for iCode = 1:length(indicatorCode)

    url = [beginUrl, char(string(indicatorCode{iCode})), '/dados?formato=csv'];

%     if ~isempty(startDate) && ~isempty(endDate)
%         url = [url, '&dataInicial=', startDate_str, '&dataFinal=' endDate_str];
%     end

    opts = delimitedTextImportOptions( ...
        'Delimiter', ';', ...
        'DataLines', 2, ...
        'VariableNames', {'date', 'values'}, ...
        'VariableTypes', {'char', 'char'});

    opt = weboptions('ContentReader', @(x, y)readtable(x, opts));
    
    dataset = webread(url, opt);

    dates   = datetime(dataset.date, 'InputFormat', 'dd/MM/yyyy');
    values  = str2double(replace(dataset.values, ',', '.'));

    ts = timetable(dates, values);
    ts.Properties.VariableNames = genvarname(string(indicatorCode{iCode}));
    ts_dates = ts.Properties.RowTimes(find((ts.Properties.RowTimes >= startDate) .* (ts.Properties.RowTimes <= endDate)));
    ts = ts(ts_dates, :);
    
    ts_final = synchronize(ts_final, ts, 'Union');
end

end