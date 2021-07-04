%% PAPER: An Unpleasant Coincidence for Monetary Policy: Risky Assets and Fiscal Limits
%
% Author: Eduardo G. C. Amaral
% Version date: May 25th 2021
%
% This routine plots the government expenses/revenues and estimates the
% fiscal rule

%% Import data from text file
clear all;
clc;

paper = 'Paper_2'; % 'Paper_2'
isFinal = false;

if strcmp(paper, 'Paper_2')
    if isFinal
        imagesFolder = '/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 2/Images/Estimation/';
    else
        imagesFolder = '/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) An Unpleasant Coincidence for Monetary Policy - Risky Assets and Fiscal Limits/Images/Calibration/';
    end
end

% Filter dates
if strcmp(paper, 'Paper_2')
    sampleYearStart = 2006;
    sampleYearEnd   = 2018;
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 59);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["CountryName", "CountryCode", "ClassificationName", "ClassificationCode", "SectorName", "SectorCode", "UnitName", "UnitCode", "Attribute", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56", "IndicatorCode", "GlobalDSDTimeSeriesCode", "VarName59"];
opts.VariableTypes = ["string", "double", "string", "double", "string", "double", "string", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string"];
opts = setvaropts(opts, [58, 59], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [4, 6, 8, 57], "TrimNonNumeric", true);
opts = setvaropts(opts, [4, 6, 8, 57], "ThousandsSeparator", ",");
opts = setvaropts(opts, [1, 3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
dataset = readtable("GFSMAB_06-13-2020 16-01-20-13_timeSeries.csv", opts);


%% Clear temporary variables
clear opts

%% Select data

% Brazil country code = 223
dsBrazil = dataset(dataset.CountryCode == 223,:);
% General government sector code = 13
dsBrazil = dsBrazil(dsBrazil.SectorCode == 13,:);
% % of GDP unit code = 1
dsBrazil = dsBrazil(dsBrazil.UnitCode == 1,:);
% Value attribute
dsBrazil = dsBrazil(ismember(dsBrazil.Attribute,"Value"),:);

firstYear   = 1972;
lastYear    = 2020;
colStart = find(ismember(dsBrazil.Properties.VariableNames, 'VarName10'));
colEnd = find(ismember(dsBrazil.Properties.VariableNames, 'VarName56'));

t = array2table(table2array(dsBrazil(:,colStart:colEnd)), ...
    'VariableNames', genvarname(string(1972:2018)) , ...
    'RowNames',string(dsBrazil.ClassificationName))

%% Sample (Expenses)

dates = (sampleYearStart:sampleYearEnd)';
colStart = find(ismember(t.Properties.VariableNames, ['x', num2str(sampleYearStart)]));
colEnd = find(ismember(t.Properties.VariableNames, ['x', num2str(sampleYearEnd)]));

itemsTotal = {'Expense', 'Compensation of employees', 'Use of goods and services', ...
         'Consumption of fixed capital', 'Other expense', 'Interest expense', ...
         'Subsidies', 'Grants expense', 'Social benefits'};

itemsGovConsumption = {'Compensation of employees', 'Use of goods and services', ...
         'Consumption of fixed capital', 'Other expense'};

itemsTransfers = {'Subsidies', 'Grants expense', 'Social benefits'};

itemsInterest = {'Interest expense'};
     
tFiltered = t(itemsTotal, colStart:colEnd)

%% Calculate Mean

tFilteredMean = tFiltered;
tFilteredMean = array2table(mean(tFilteredMean.Variables,2,'omitnan'), ...
    'VariableNames', {'Mean'} , ...
    'RowNames', itemsTotal);

disp(tFilteredMean)

datesVector = datetime(str2double(cellfun(@(x) x(2:end), tFiltered.Properties.VariableNames, 'UniformOutput', false)), 1, 1)';
govExpense = tFiltered('Expense',:).Variables';
tGovExpense = timetable(datesVector, govExpense, 'VariableNames', {'govExpense'});


%% Plot results: Expenses

yExpense    = tFiltered{"Expense", :};
yGovCons    = sum(tFiltered{itemsGovConsumption, :}, 1, 'omitnan')';
yTransfers  = sum(tFiltered{itemsTransfers, :}, 1, 'omitnan')';
yInterest   = sum(tFiltered{itemsInterest, :}, 1, 'omitnan')';

f = figure;
pub_GraphSetInterpreter('latex');
y = [yGovCons, yGovCons + yTransfers, yGovCons + yTransfers + yInterest];
p = plot(dates, y);
p(1).LineStyle = '-';
p(2).LineStyle = '--';
p(3).LineStyle = '-.';
p(1).LineWidth = 3.5;
p(2).LineWidth = 3.5;
p(3).LineWidth = 3.5;
p(1).Marker = '*';
p(2).Marker = '*';
p(3).Marker = '*';
p(1).MarkerSize = 10;
p(2).MarkerSize = 10;
p(3).MarkerSize = 10;
set(gca,'FontSize', 14);
ylabel('\% of GDP');
xlabel('date');
xlim([sampleYearStart sampleYearEnd]);
ylim([0 ceil(1.1*max(y(:)) / 5) * 5]);
grid('on');
set(gca,'LooseInset',get(gca,'TightInset'));
title('Stacked Brazilian government expenses');
[h, icons] = legend('Goverment consumption', 'Transfer expenses', 'Interest expenses', 'location', 'best');
% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',0.0001);

simGraphName = 'actual_governmentExpenses.png';
set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);


%% Sample (Revenue)

dates = (sampleYearStart:sampleYearEnd)';
colStart = find(ismember(t.Properties.VariableNames, ['x', num2str(sampleYearStart)]));
colEnd = find(ismember(t.Properties.VariableNames, ['x', num2str(sampleYearEnd)]));

itemsTotal = {'Revenue', 'Tax revenue', 'Social contributions', ...
         'Grants revenue', 'Other revenue'};

itemsTaxes = {'Tax revenue'};

itemsOtherRevenues = {'Social contributions', ...
         'Grants revenue', 'Other revenue'};

     
tFiltered = t(itemsTotal, colStart:colEnd)

%% Calculate Mean

tFilteredMean = tFiltered;
tFilteredMean = array2table(mean(tFilteredMean.Variables,2,'omitnan'), ...
    'VariableNames', {'Mean'} , ...
    'RowNames', itemsTotal);

disp(tFilteredMean)

datesVector = datetime(str2double(cellfun(@(x) x(2:end), tFiltered.Properties.VariableNames, 'UniformOutput', false)), 1, 1)';
govRevenue = tFiltered('Revenue',:).Variables';
tGovRevenue = timetable(datesVector, govRevenue, 'VariableNames', {'govRevenue'});

%% Plot results: Revenues

yRevenue    = tFiltered{"Revenue", :};
yTaxes    = sum(tFiltered{itemsTaxes, :}, 1, 'omitnan')';
yOtherRevenues  = sum(tFiltered{itemsOtherRevenues, :}, 1, 'omitnan')';

f = figure;
pub_GraphSetInterpreter('latex');
y = [yTaxes, yTaxes + yOtherRevenues];
p = plot(dates, y);
p(1).LineStyle = '-';
p(2).LineStyle = '--';
p(1).LineWidth = 3.5;
p(2).LineWidth = 3.5;
p(1).Marker = '*';
p(2).Marker = '*';
p(1).MarkerSize = 10;
p(2).MarkerSize = 10;
set(gca,'FontSize', 14);
ylabel('\% of GDP');
xlabel('date');
xlim([sampleYearStart sampleYearEnd]);
ylim([0 ceil(1.1*max(y(:)) / 5) * 5]);
grid('on');
set(gca,'LooseInset',get(gca,'TightInset'));
title('Stacked Brazilian government revenues');
[h, icons] = legend('Tax revenues', 'Other revenues', 'location', 'best');
% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',0.0001);

simGraphName = 'actual_governmentRevenues.png';
set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%% Estimate gammaTau

% Brazil country code = 223
dsBrazil = dataset(dataset.CountryCode == 223,:);
% General government sector code = 13
dsBrazil = dsBrazil(dsBrazil.SectorCode == 13,:);
% % of GDP unit code = 1 // domestic currency = NaN
dsBrazil = dsBrazil(dsBrazil.UnitCode == 1,:);
% Value attribute
dsBrazil = dsBrazil(ismember(dsBrazil.Attribute,"Value"),:);

firstYear   = 1972;
lastYear    = 2020;
colStart = find(ismember(dsBrazil.Properties.VariableNames, 'VarName10'));
colEnd = find(ismember(dsBrazil.Properties.VariableNames, 'VarName56'));

t = array2table(table2array(dsBrazil(:,colStart:colEnd)), ...
    'VariableNames', genvarname(string(1972:2018)) , ...
    'RowNames',string(dsBrazil.ClassificationName));

% Sample (Expenses)
dates = (sampleYearStart:sampleYearEnd)';
colStart = find(ismember(t.Properties.VariableNames, ['x', num2str(sampleYearStart)]));
colEnd = find(ismember(t.Properties.VariableNames, ['x', num2str(sampleYearEnd)]));

itemsTotal = {'Expense', 'Revenue'};
tFiltered = t(itemsTotal, colStart:colEnd);
tFiltered = array2timetable( ...
    tFiltered.Variables', ...
    'RowTimes', datetime(str2double(cellfun(@(x) x(2:end), tFiltered.Properties.VariableNames, 'UniformOutput', false)), 1,1)', ...
    'VariableNames', tFiltered.Properties.RowNames);

% Get Debt to GDP
% 13762	Dívida bruta do governo geral - Metodologia utilizada a partir de 2008
% 1207	PIB em R$ correntes
% 22099	PIB a preços de mercado SCN-2010
seriesCode = 4537;
tDebtToGDP = pub_Database_SGS_Get(seriesCode, datetime(sampleYearStart,1,1), datetime(sampleYearEnd,1,1));
tDebtToGDP.Properties.VariableNames = {'debtToGDP'};

seriesCode = 1207;
tNominalGDP = pub_Database_SGS_Get(seriesCode, datetime(sampleYearStart,1,1), datetime(sampleYearEnd,1,1));
tNominalGDP.Properties.VariableNames = {'nominalGDP'};

seriesCode = 22099;
tRealGDP = pub_Database_SGS_Get(seriesCode, datetime(sampleYearStart,1,1), datetime(sampleYearEnd,1,1));
tRealGDP.Properties.VariableNames = {'realGDP'};

% Consolidate
tReg = synchronize(tFiltered, tDebtToGDP, 'yearly', 'mean');
tReg = synchronize(tReg, tNominalGDP, 'yearly', 'mean');
tReg = synchronize(tReg, tRealGDP, 'yearly', 'mean');
tReg.debtToGDP_LAG1     = lagmatrix(tReg{:, 'debtToGDP'}, 1);
tReg.nominalGDP_LAG1    = lagmatrix(tReg{:, 'nominalGDP'}, 1);
tReg.realGDP_LAG1       = lagmatrix(tReg{:, 'realGDP'}, 1);
tReg.revenue_LAG1       = lagmatrix(tReg{:, 'Revenue'}, 1);
tReg.expense_LAG1       = lagmatrix(tReg{:, 'Expense'}, 1);

% % Fit Model 1
% disp('Regression: Revenue ~ debtToGDP(-1)');
% y = log(tReg{:, 'Revenue'} ./ mean(tReg{:, 'Revenue'}));
% X = log(tReg{:, 'debtToGDP_LAG1'} ./ mean(tReg{:, 'debtToGDP'}));
% mdl = fitlm(X, y)
% plot(mdl)

% Fit Model 2 (This is the specification of the paper)
disp('Regression: Revenue ~ Revenue(-1) + debtToGDP(-1)');
y = log(tReg{:, 'Revenue'} ./ mean(tReg{:, 'Revenue'}));
X = [log(tReg{:, 'revenue_LAG1'} ./ mean(tReg{:, 'Revenue'})), log(tReg{:, 'debtToGDP_LAG1'} ./ mean(tReg{:, 'debtToGDP'}))];
mdl = fitlm(X, y, 'RobustOpts','off')
f = figure;
plot(mdl);
f = figure;
plotResiduals(mdl, 'probability');

% % Fit Model 3
% disp('Regression: Expense ~ Expense(-1)');
% y = log(tReg{:, 'Expense'} ./ mean(tReg{:, 'Expense'}));
% X = [log(tReg{:, 'expense_LAG1'} ./ mean(tReg{:, 'Expense'}))];
% mdl = fitlm(X, y)
% plot(mdl)

% % Fit Model 4
% disp('Regression: Expense ~ Expense(-1) + Output(-1)');
% y = log(tReg{:, 'Expense'} ./ mean(tReg{:, 'Expense'}));
% X = [log(tReg{:, 'expense_LAG1'} ./ mean(tReg{:, 'Expense'})), log(tReg{:, 'realGDP_LAG1'} ./ mean(tReg{:, 'realGDP'}))];
% mdl = fitlm(X, y)
% plot(mdl)
