%% Fetch data from IMF World Economic Outlook

function result = pub_Database_WEO_Fetch(WEODatabase, WEOSubjectCode, ISO, returnedCols)

%%%%%% DESCRIPTION
% WEODatabase       = WEO table obtained from pub_Database_WEO_Load
% WEOSubjectCode    = code of the indicator
% ISO               = country 3-letter ISO code
% returnedCols      = names of the columns to be returned
%%%%%%

% Search indicator
idx1 = string(WEODatabase.WEOSubjectCode) == WEOSubjectCode;

% Search country
if isempty(ISO)
    idx2 = ones(size(WEODatabase.ISO,1),1);
else
    idx2 = string(WEODatabase.ISO) == ISO;
end

% Combine filters
idx = repmat(idx1,1,size(idx2,2)) .* idx2;

idxWEO      = ismember(WEODatabase.Properties.VariableNames, genvarname(string(returnedCols)));
idxReturn   = ismember(genvarname(string(returnedCols)), WEODatabase.Properties.VariableNames(idxWEO));

idxISO          = ismember(WEODatabase.Properties.VariableNames, 'ISO');
idxSubjectCode  = ismember(WEODatabase.Properties.VariableNames, 'WEOSubjectCode');

% Final result
dates = datetime(returnedCols(idxReturn),1,1)';
[tbIdx_row, ~] = find(idx);
[~, tbIdxWEO_col] = find(idxWEO);

if length(ISO) == 1
    seriesName = WEODatabase{tbIdx_row, idxSubjectCode};
else
    seriesName = WEODatabase{tbIdx_row, idxISO};
end

tb = array2table(cell2mat(WEODatabase{tbIdx_row, tbIdxWEO_col})', 'VariableNames', seriesName);

if isempty(tb)
    disp('No results');
    result = [];
    return
end

result = table2timetable(tb, 'RowTimes', dates);

end