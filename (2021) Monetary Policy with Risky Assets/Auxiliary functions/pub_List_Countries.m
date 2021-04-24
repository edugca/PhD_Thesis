%% Create table with countries information
function tb_countries = pub_List_Countries()

% Load countries lists
t_ISOCountryCodes = pub_Database_WEO_CountryCodes_Load();
WEODatabase = pub_Database_WEO_Load();

advancedCountries = pub_Database_WEO_CountryGroups_Load('Advanced Economies');
emergingCountries = pub_Database_WEO_CountryGroups_Load('Emerging and Developing Economies');

% Adjust names of the countries
countryCorresp = unique(string([WEODatabase.Country WEODatabase.ISO]), 'rows');
emergingCountries{string(emergingCountries) == 'Iran'} = 'Islamic Republic of Iran';
emergingCountries(string(emergingCountries) == 'Kosovo') = []; % There is no ISO code for Kosovo
emergingCountries = countryCorresp(cellfun(@(x) find(x == countryCorresp(:,1)), emergingCountries), 2);
advancedCountries = countryCorresp(cellfun(@(x) find(x == countryCorresp(:,1)), advancedCountries), 2);

% Select countries
countryList = [advancedCountries; emergingCountries];
countryStatus = [repmat("Advanced", length(advancedCountries), 1); repmat("Emerging", length(emergingCountries), 1)];

% Make lists of countries and their codes
countryCodes_2L = t_ISOCountryCodes(ismember(t_ISOCountryCodes.ISO3, countryList), 2).Variables;
countryCodes_3L = t_ISOCountryCodes(ismember(t_ISOCountryCodes.ISO3, countryList), 3).Variables;
countryNames    = t_ISOCountryCodes(ismember(t_ISOCountryCodes.ISO3, countryList), 1).Variables;
countryStatus   = string();
countryStatus(ismember(countryCodes_3L, emergingCountries)) = 'Emerging';
countryStatus(ismember(countryCodes_3L, advancedCountries)) = 'Advanced';
countryStatus = countryStatus';
tb_countries = table(countryNames, countryCodes_2L, countryCodes_3L, countryStatus, ...
    'VariableNames', {'Name', 'ISO2', 'ISO3', 'Status'}, 'RowNames', countryCodes_3L);

% Add list of alternative names
tb_countries.Name_Alternative_1 = strings(size(tb_countries,1), 1);
tb_countries{'GBR','Name_Alternative_1'} = "United Kingdom";
tb_countries{'IRN','Name_Alternative_1'} = "Iran";
tb_countries{'USA','Name_Alternative_1'} = "United States";

end