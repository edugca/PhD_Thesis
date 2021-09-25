%% Fetch data of CDS from Bloomberg dataset

function result = pub_Database_LCCS_Fetch()

%%%%%% DESCRIPTION
% sourceName       = 'Bloomberg'
% indicatorName    = name of the variable to be exported
%%%%%%

pathData = pub_Path('Data/LCCS');
filename = 'cip_dataset_v3.csv';
filePath = [pathData filesep filename];

% Specify sheet and range
opts = detectImportOptions(filePath);
opts.VariableTypes(1:2) = {'string'};
opts.VariableTypes(3) = {'string'};
opts.VariableTypes(4) = {'string'};
opts.VariableTypes(5:end) = {'double'};

% Import the data
dataset = readtable(filePath, opts);

% Format variables
dataset.date = datetime(datevec(dataset.date),'Format','ddMMMyyyy');

ctrTable = array2table([    "AUD", "AU";
                "BRL", "BR";
                "CAD", "CA";
                "CHF", "CH";
                "CLP", "CL";
                "CNY", "CN";
                "COP", "CO";
                "DKK", "DK";
                "EUR", "EU";
                "GBP", "GB";
                "HUF", "HU";
                "IDR", "ID";
                "ILS", "IL";
                "INR", "IN";
                "JPY", "JP";
                "KRW", "KR";
                "MXN", "MX";
                "MYR", "MY";
                "NOK", "NO";
                "NZD", "NZ";
                "PEN", "PE";
                "PHP", "PH";
                "PLN", "PL";
                "RUB", "RU";
                "SEK", "SE";
                "THB", "TH";
                "TRY", "TR";
                "ZAR", "ZA";
           ]);
ctrTable.Properties.VariableNames = {'currency', 'country'};
curTable.Properties.VariableNames = {'currency'};
dataset = innerjoin(dataset, ctrTable);

dataset.group = [];
dataset.diff_y = [];
dataset.rho = [];
dataset = dataset(strcmp(dataset.tenor, "5y"), :);
dataset.tenor = [];
dataset.currency = [];

wideDataset = unstack(dataset, 'cip_govt', 'country');
wideDataset = table2timetable(wideDataset);

% Final result
result = wideDataset;

end