%% Fetch data from IMF World Economic Outlook

function result = pub_Database_WEO_CountryCodes_Load()

pathData = pub_Path('/Users/Eduardo/OneDrive/MATLAB/My Library/Database/Data/WEO', 'C:\');

filename = 'Country Codes.xlsx';
filePath = [pathData filesep filename];


% Specify sheet and range
opts = detectImportOptions(filePath);
opts = setvartype(opts, 'char');

% Import the data
listData = readtable(filePath, opts);

% Final result
result = listData;

end