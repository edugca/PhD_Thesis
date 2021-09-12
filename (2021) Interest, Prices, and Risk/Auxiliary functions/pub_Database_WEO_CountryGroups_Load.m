%% Fetch data from IMF World Economic Outlook

function result = pub_Database_WEO_CountryGroups_Load(groupName)

pathData = pub_Path('/Users/Eduardo/OneDrive/MATLAB/My Library/Database/Data/WEO', 'C:\');
filename = 'Country Classification.xlsx';
filePath = [pathData filesep filename];

% Specify sheet and range
opts = detectImportOptions(filePath);
opts = setvartype(opts, 'char');

% Import the data
listData = readtable(filePath, opts);

% Select group
groupCountries = listData.(genvarname(groupName));

% Remove empty entries
groupCountries = groupCountries(cellfun(@(x) ~isempty(x), groupCountries));

% Final result
result = groupCountries;

end