%% Fetch data from IMF World Economic Outlook

function result = pub_Database_WEO_Load(varargin)

%%%%%% DESCRIPTION
% Year    = i.e '2010'
% Month   = i.e. 'Apr', 'Oct' or 'Sep'
%%%%%%

% Latest WEO
weo_year = '2020';
weo_month = 'Oct';

if length(varargin) > 1
    for i = 1:length(varargin)
        if strcmp(varargin{i}, 'Year')
            weo_year = num2str(varargin{i+1});
        end
        if strcmp(varargin{i}, 'Month')
            weo_month = num2str(varargin{i+1});
        end
    end
end

pathData = pub_Path(['/Users/Eduardo/OneDrive/MATLAB/My Library/Database/Data/WEO/' weo_year], 'C:\');

filename = ['WEO' weo_month weo_year 'all.xls'];
filePath = [pathData filesep filename];

% Specify sheet and range
if ~isfile(filePath)
    disp('WEO database file could not be found');
    result = [];
    return 
end
opts = detectImportOptions(filePath, 'FileType', 'delimitedText');
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
opts = setvartype(opts, 'char');

% Import the data
surveyData = readtable(filePath, opts);

% Convert to numbers
colsNumbers = startsWith(surveyData.Properties.VariableNames, 'x');
surveyData{:,colsNumbers} = cellfun(@(x) str2double(x), table2cell(surveyData(:,colsNumbers)), 'UniformOutput',false);

% Final result
result = surveyData;

end