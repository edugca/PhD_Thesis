%% Fetch data of CDS from Bloomberg dataset

function result = pub_Database_CDS_Fetch(sourceName)

%%%%%% DESCRIPTION
% sourceName       = 'Bloomberg'
% indicatorName    = name of the variable to be exported
%%%%%%

pathData = pub_Path('Data/CDS');
filename = 'CDS.xlsx'; %'CDS.xlsx'
filePath = [pathData filesep filename];

% Specify sheet and range
opts = detectImportOptions(filePath, 'Sheet', sourceName, 'NumHeaderLines', 6);
opts.VariableTypes(2:end) = {'double'};
opts.VariableNamesRange = 'A4';

% Import the data
dataset = readtimetable(filePath, opts, "UseExcel", false, 'ReadVariableNames', true);

% Final result
result = dataset;

end