%%%%%%%%% Amaral(2021): Monetary Policy with Risky Assets

% This script plots figures 1, 21, and 22

%% Set up

clear all; clc;

isFinal = false;

if isFinal
    pathData                = pub_Path('/Users/Eduardo/OneDrive/MATLAB/My Projects/Natural Rate of Interest/Practitioner Methods', 'C:\');
    pathSaved               = pub_Path('/Users/Eduardo/OneDrive/MATLAB/My Library/Database/Data/My Projects/', 'C:\');
    
    imagesFolder            = pub_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Images/Filters/', 'C:\');
    tablesFolder            = pub_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Tables/Filters/', 'C:\');
else
    pathData                = pub_Path('Data/');
    pathSaved               = pub_Path('Filters/Saved/');
    
    imagesFolder            = pub_Path('Filters/Images/');
    tablesFolder            = pub_Path('Filters/Tables/');
    
end

updateFRED      = false;
updateSPF       = false;
updateWEO       = false;
updateIMF       = false;
updateFOCUS     = false;
updateOutputGap = false;

originalData = timetable();

%% OECD countries + Others

load([pathData, 'isoCountryCodes.mat']);

countryCodesOECD_2L = ['AT'; 'AU'; 'BE'; 'CA'; 'CZ'; 'DK'; 'EE'; 'FI'; 'FR'; 'DE'; 'GR'; 'HU'; 'IS'; 'IE'; 'IL'; 'IT'; 'JP'; 'KR'; 'LU'; 'MX'; 'NL'; 'NZ'; 'NO'; 'PL'; 'PT'; 'SK'; 'SI'; 'ES'; 'SE'; 'CH'; 'TR'; 'GB'; 'US'; 'EZ'; 'EU'; 'O1'; 'G7'];
countryCodesALL_2L     = [countryCodesOECD_2L; 'BR'; 'CL'; 'CN'];

countryCodes_2L = t_ISOCountryCodes(ismember(t_ISOCountryCodes.Alpha2code, countryCodesALL_2L), 2).Variables;
countryCodes_3L = t_ISOCountryCodes(ismember(t_ISOCountryCodes.Alpha2code, countryCodesALL_2L), 3).Variables;
countryNames = t_ISOCountryCodes(ismember(t_ISOCountryCodes.Alpha2code, countryCodesALL_2L), 1).Variables;

%% Download series from IMF-WEO

% Load countries lists
t_ISOCountryCodes = pub_Database_WEO_CountryCodes_Load();
WEODatabase = pub_Database_WEO_Load();

% Load my countries' list
tb_countries        = pub_List_Countries();
countryNames        = tb_countries.Name(:);
countryCodes_2L     = tb_countries.ISO2(:);
countryCodes_3L     = tb_countries.ISO3(:);
countryList         = tb_countries.ISO3(:);
countryStatus       = tb_countries.Status(:);
advancedCountries   = tb_countries{find(strcmp(tb_countries{:,'Status'}, 'Advanced')),'ISO3'};
emergingCountries   = tb_countries{find(strcmp(tb_countries{:,'Status'}, 'Emerging')),'ISO3'};

ii = 0;
series_names        = {};
series_codes        = {};

startDate = 1900;
endDate   = 2100;

% Nominal GDP in US Dollars (NGDPD)
ii = ii + 1;
series_names{ii}        = 'NY_';
series_codes{ii}        = 'NGDPD';

% Gross domestic product, constant prices (NGDP_R)
ii = ii + 1;
series_names{ii}        = 'RY_';
series_codes{ii}        = 'NGDP_R';

% Employment, millions of people (LE)
ii = ii + 1;
series_names{ii}        = 'EMP_';
series_codes{ii}        = 'LE';

% Fetch series
if updateWEO
    originalData_WEO = timetable();
    for ii = 1:length(series_names)
        for jj = 1:length(countryCodes_3L)
            aux = pub_Database_WEO_Fetch(...
                WEODatabase, series_codes{ii}, countryCodes_3L{jj}, ...
                startDate:endDate);

            % Check whether it was successful
            if istimetable(aux)
                aux.Properties.VariableNames = {[series_names{ii} countryCodes_2L{jj}]};
                originalData_WEO = synchronize(originalData_WEO, aux);
            end
        end
    end
    
    save([pathData, 'WEO.mat'], 'originalData_WEO');
else
    load([pathData, 'WEO.mat']);
end

originalData = synchronize(originalData, originalData_WEO);

% Filter countries: Find the cutN largest economies at date cutDate
cutN = 20;
cutDate = datetime(2019,1,1);

listAdvanced = tb_countries(ismember(tb_countries.Status,{'Advanced'}),'ISO2');
listEmerging = tb_countries(ismember(tb_countries.Status,{'Emerging'}),'ISO2');

%%%% Advanced economies
cutVars = join([repmat('NY_',size(listAdvanced,1), 1) string(listAdvanced.Variables)], '');
idxVars = ismember(originalData_WEO.Properties.VariableNames, cutVars);
advancedData = originalData_WEO(cutDate,idxVars);
[~, largestOrder] = sort(advancedData(cutDate,:).Variables, ...
    'descend', 'MissingPlacement', 'last');
sortedValues = advancedData(cutDate,largestOrder);
sortedValues = sortedValues(:, 1:cutN);

selAdvancedCountries_2L = cellfun(@(x) x(end-1:end), ...
    sortedValues.Properties.VariableNames, 'UniformOutput', false)';
[~, idx] = ismember(t_ISOCountryCodes.ISO2, selAdvancedCountries_2L);
tb_countries.Properties.RowNames = tb_countries.ISO2;
selAdvancedCountries_3L = tb_countries(selAdvancedCountries_2L, 3).Variables;
selAdvancedCountries_Names = tb_countries(selAdvancedCountries_2L, 1).Variables;
countryWeight_Advanced = table(sortedValues.Variables', 'VariableNames', {'Weight'}, 'RowNames', selAdvancedCountries_2L);
%%%%%%%%%%%%%%%%%%%%%%

%%%% Emerging economies
cutVars = join([repmat('NY_',size(listEmerging,1), 1) string(listEmerging.Variables)], '');
idxVars = ismember(originalData_WEO.Properties.VariableNames, cutVars);
emergingData = originalData_WEO(cutDate,idxVars);
[~, largestOrder] = sort(emergingData(cutDate,:).Variables, ...
    'descend', 'MissingPlacement', 'last');
sortedValues = emergingData(cutDate,largestOrder);
sortedValues = sortedValues(:, 1:cutN);

selEmergingCountries_2L = cellfun(@(x) x(end-1:end), ...
    sortedValues.Properties.VariableNames, 'UniformOutput', false)';
[~, idx] = ismember(t_ISOCountryCodes.ISO2, selEmergingCountries_2L);
tb_countries.Properties.RowNames = tb_countries.ISO2;
selEmergingCountries_3L = tb_countries(selEmergingCountries_2L, 3).Variables;
selEmergingCountries_Names = tb_countries(selEmergingCountries_2L, 1).Variables;
countryWeight_Emerging = table(sortedValues.Variables', 'VariableNames', {'Weight'}, 'RowNames', selEmergingCountries_2L);
%%%%%%%%%%%%%%%%%%%%%%

%%%% All economies
selAllCountries_2L = [selAdvancedCountries_2L; selEmergingCountries_2L];
selAllCountries_3L = [selAdvancedCountries_3L; selEmergingCountries_3L];
selAllCountries_Names = [selAdvancedCountries_Names; selEmergingCountries_Names];
countryWeight_All = [countryWeight_Advanced; countryWeight_Emerging];
%%%%%%%%%%%%%%%%%%%%%%

%% Download series from IMF-IFS

ii = 0;
series_names        = {};
series_databases    = {};
series_codes        = {};
series_frequencies  = {};

startDate = '1900';
endDate   = '2100';

% Central Bank Policy Rate (FPOLM_PA)
ii = ii + 1;
series_names{ii}        = 'NR_';
series_databases{ii}    = 'IFS';
series_codes{ii}        = 'FPOLM_PA';
series_frequencies{ii}  = 'M';

% Money Market Rate (FIMM_PA)
ii = ii + 1;
series_names{ii}        = 'MM_';
series_databases{ii}    = 'IFS';
series_codes{ii}        = 'FIMM_PA';
series_frequencies{ii}  = 'M';

% Treasury Bill Rate (FITB_PA)
ii = ii + 1;
series_names{ii}        = 'TB_';
series_databases{ii}    = 'IFS';
series_codes{ii}        = 'FITB_PA';
series_frequencies{ii}  = 'M';
 
% Prices, Consumer Price Index, All items, Index (PCPI_IX)
ii = ii + 1;
series_names{ii}        = 'CPI_';
series_databases{ii}    = 'IFS';
series_codes{ii}        = 'PCPI_IX';
series_frequencies{ii}  = 'M';

% National Accounts, Expenditure, Gross Domestic Product, Real, Domestic Currency (NGDP_R_XDC)
ii = ii + 1;
series_names{ii}        = 'Y_';
series_databases{ii}    = 'IFS';
series_codes{ii}        = 'NGDP_R_XDC';
series_frequencies{ii}  = 'Q';

% Fetch series
if updateIMF
    originalData_IMF = timetable();
    for ii = 1:length(series_names)
        for jj = 1:length(selAllCountries_2L)
            aux = pub_Database_IMF_Fetch( ...
                series_databases{ii}, series_codes{ii}, ...
                selAllCountries_2L{jj}, series_frequencies{ii}, ...
                startDate, endDate);

            % Check whether it was successful
            if istimetable(aux)
                aux.Properties.VariableNames = {[series_names{ii} selAllCountries_2L{jj}]};
                originalData_IMF = synchronize(originalData_IMF, aux);
            end
        end
    end
    
    save([pathData, 'IFS.mat'], 'originalData_IMF');
else
    load([pathData, 'IFS.mat']);
end

originalData = synchronize(originalData, originalData_IMF);

%% Treat series

% Quarterly the dataset
ts = retime(originalData,'quarterly','mean');

%%%%%%%%%% OECD and others

% CPI: Quarterly annualized % change
selectedVars = startsWith(ts.Properties.VariableNames, 'CPI_');
for jj=1:length(ts.Properties.VariableNames)
    if selectedVars(jj) == 1
        ts.(ts.Properties.VariableNames{jj}) = 100 * [nan; (ts.(ts.Properties.VariableNames{jj})(2:end) ./ ts.(ts.Properties.VariableNames{jj})(1:end-1)).^4 - 1];
    end
end

% CPI Forecast: Quarterly annualized % change
nHorizon = 4;
pMax = 4;
dMax = 0;
qMax = 0;
seasonality = 0;
sarMax = {};
smaMax = {};
nEstObs = 80; % number of last periods to use in the estimation
selectedVars = any(  [startsWith(ts.Properties.VariableNames, 'CPI_'); ...
                    startsWith(ts.Properties.VariableNames, 'Y_'); ...
                    startsWith(ts.Properties.VariableNames, 'IPEA_'); ...
                    startsWith(ts.Properties.VariableNames, 'IFI_')]);
iSpec = 0;
for jj=1:length(ts.Properties.VariableNames)
    if selectedVars(jj) == 1
        
        % Display series
        disp(ts.Properties.VariableNames{jj});
        
        yMdl = ts.(ts.Properties.VariableNames{jj});
        nObs = sum(~isnan(yMdl));
        
        if nObs > 1
            bestMdl = pub_ARIMABestFit(yMdl(end-nEstObs+1:end), ...
                pMax, dMax, qMax, 'bic', seasonality, sarMax, smaMax);

            mdl = arima(bestMdl(1), bestMdl(2), bestMdl(3));
            mdl.Seasonality = seasonality;
            estMdl = estimate(mdl, yMdl(end-nEstObs+1:end));

            yForecast = pub_ModelRollingForecast(estMdl, nHorizon, yMdl);
            yForecast = real( (prod((1+yForecast./100).^(1/4), 2) - 1).*100 );
        else
           bestMdl = 'none';
           yForecast = NaN(length(yMdl),1);
        end
        
        iSpec = iSpec + 1;
        varName = [ts.Properties.VariableNames{jj}, '_Forecast_1Y'];
        varClassName = extractBefore(ts.Properties.VariableNames{jj},'_');
        eval(['mdl_', varClassName, '_Forecast_Specs{iSpec,1} = varName;']);
        eval(['mdl_', varClassName, '_Forecast_Specs{iSpec,2} = bestMdl;']);
        yTimeTable = timetable(ts.Properties.RowTimes, yForecast, 'VariableNames', {varName});
        
        ts = synchronize(ts, yTimeTable);
    end
end
mdl_CPI_Forecast_Specs = mdl_CPI_Forecast_Specs(cellfun(@(x) ~isempty(x), mdl_CPI_Forecast_Specs(:,1)),:);
mdl_Y_Forecast_Specs = mdl_Y_Forecast_Specs(cellfun(@(x) ~isempty(x), mdl_Y_Forecast_Specs(:,1)),:);

% Output: take the log
selectedVars = startsWith(ts.Properties.VariableNames, 'Y_');
for jj=1:length(ts.Properties.VariableNames)
    if selectedVars(jj) == 1
    ts.(['log', ts.Properties.VariableNames{jj}]) = log(ts.(ts.Properties.VariableNames{jj}));
    ts.(['pctChg_', ts.Properties.VariableNames{jj}])   = [nan; ( ( ts.(ts.Properties.VariableNames{jj})(2:end) ./ ts.(ts.Properties.VariableNames{jj})(1:end-1) ).^4 - 1 ) * 100];
    ts.(['YoY_' ts.Properties.VariableNames{jj}])     = (movavg(ts.(ts.Properties.VariableNames{jj}),'linear',4) ./ lagmatrix(movavg(ts.(ts.Properties.VariableNames{jj}),'linear',4),4) - 1) * 100;
    end
end

%%%%%%%%%%

% Export data
fileDataPath = [pathData, 'rStar_dataset.xlsx'];
writetimetable(ts, fileDataPath);

% Auxiliary variables
nPeriods = size(ts,1);

%% SAVE/LOAD DATA

%save([pathSaved, 'rStar_workspace.mat'])

load([pathSaved, 'rStar_workspace.mat']);


%% 1) Filters: HP, Baxter-King, Christiano-Fitzgerald (All Countries)

close all; clc;

ts_All_Filters = struct();

nCountries  = sum(startsWith(ts.Properties.VariableNames, 'NR_'));
iCountry    = 0;

% Plot individual countries
fAdvanced = figure;
fEmerging = figure;

previousInterpreter = pub_GraphSetInterpreter('latex');

% Set layout and reduce empty spaces in the graph
tlAdvanced = tiledlayout(fAdvanced, 'flow');
tlAdvanced.TileSpacing  = 'compact';
tlAdvanced.Padding      = 'compact';
legendInFirstSubPlot = false;
if legendInFirstSubPlot
   legend_axes_advanced = nexttile(tlAdvanced); 
end

tlEmerging = tiledlayout(fEmerging, 'flow');
tlEmerging.TileSpacing  = 'compact';
tlEmerging.Padding      = 'compact';
legendInFirstSubPlot = false;
if legendInFirstSubPlot
   legend_axes_emerging = nexttile(tlEmerging); 
end

tDataSample = timetable(ts.Properties.RowTimes);

selForPloting = {};
for jj=1:length(selAllCountries_Names)
    % Check if there is data available
    if ismember(['TB_', selAllCountries_2L{jj} ], ts.Properties.VariableNames) ...
         && ismember(['CPI_', selAllCountries_2L{jj} '_Forecast_1Y'], ts.Properties.VariableNames) ...
     
        % Check whether there is a valid CPI forecasting model
        idxMdl = find(['CPI_', selAllCountries_2L{jj} '_Forecast_1Y'] == string(mdl_CPI_Forecast_Specs(:,1)));
        if sum(mdl_CPI_Forecast_Specs{idxMdl,2})> 0
     
            % Names of specifications
            specNames = {   selAllCountries_Names{jj} };
            if length(specNames{:}) > 30
                nameAlt = tb_countries{find(strcmp(tb_countries.Name, specNames)), 'Name_Alternative_1'};
                if length(nameAlt) <  length(specNames{:})
                    specNames = nameAlt;
                end
            end
            selForPloting = [selForPloting; selAllCountries_2L{jj}];
            
                        
            % Nominal rate
            nrList       = [    ts.(['TB_', selAllCountries_2L{jj}])
                            ];

            % Inflation expectation rate
            piiLRList    = [ lagmatrix(ts.(['CPI_', selAllCountries_2L{jj}, '_Forecast_1Y']),-4)
                            ];

            % Check sample
            tDataSample.([selAllCountries_2L{jj} '_NR'])    = ts.(['TB_', selAllCountries_2L{jj}]) ;
            tDataSample.([selAllCountries_2L{jj} '_CPI'])   = ts.(['CPI_', selAllCountries_2L{jj}]) ;
            
            % Apply filters
            config_plotAllFilters       = false;
            config_plotFiltersRange     = false;
            config_plotFiltersRangeAllCountries = true;
            iCountry = iCountry + 1;
            method_Filters;
            
            currentStatus = tb_countries{find(selAllCountries_2L{jj} == string(tb_countries{:,{'ISO2'}})), 'Status'};
            plot_Filters_AllCountries;
            
            ts_All_Filters.(selAllCountries_2L{jj}) = ts_Filters;
        end
    end
end

% Custom font
set(fAdvanced.Children.Children, 'FontSize', 12);
set(fEmerging.Children.Children, 'FontSize', 12);
set(fAdvanced.Children.Children, 'FontWeight', 'bold');
set(fEmerging.Children.Children, 'FontWeight', 'bold');

% Plot or Resize legend
if legendInFirstSubPlot
    legend_axes_advanced.Legend.FontSize = 14;
    legend_axes_emerging.Legend.FontSize = 14;
else
    leg = legend(tlAdvanced.Children(1), 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
    leg.FontSize = 14;
    
    leg = legend(tlEmerging.Children(1), 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';
    leg.FontSize = 14;
end

% Plot zero horizontal axis
advancedAxes = findall(fAdvanced,'type','axes');
emergingAxes = findall(fEmerging,'type','axes');
for iGraph = 1:(length(advancedAxes)-legendInFirstSubPlot)
   pub_GraphDrawZeroAxis(advancedAxes(iGraph)); 
end
for iGraph = 1:(length(emergingAxes)-legendInFirstSubPlot)
   pub_GraphDrawZeroAxis(emergingAxes(iGraph)); 
end

% Synchronize X axis
%graphDates = [datetime(2000,1,1) datetime(2020,1,1)];
%graphXTicks = graphDates(1):calyears(5):graphDates(end);
%edu_GraphSynchronizeAxes(fAdvanced, graphDates, [], graphXTicks, []);
%edu_GraphSynchronizeAxes(fEmerging, graphDates, [], graphXTicks, []);

% Save graph of individual advanced economies
f = figure(fAdvanced);
strSuptitle = ['Advanced Economies: Real neutral interest rate'];
%edu_Suptitle(strSuptitle, 'FontSize', 16);
graphName = 'rStar_Filters_Advanced.png';
set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, graphName]);

% Save graph of individual emerging economies
f = figure(fEmerging);
strSuptitle = ['Emerging Economies: Real neutral interest rate'];
%edu_Suptitle(strSuptitle, 'FontSize', 16);
graphName = 'rStar_Filters_Emerging.png';
set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, graphName]);

%%%%%%%%% Report data samples used
tSampleDes = array2table(strings(length(selForPloting), 6), ...
    'RowNames', selForPloting, ...
    'VariableNames', {'Name', 'Status', 'NR_Start', 'NR_End', 'CPI_Start', 'CPI_End'});
for iCountry = 1:length(selForPloting)
    tNR_Series = tDataSample(:, [selForPloting{iCountry}, '_NR']);
    tNR_Series = rmmissing(tNR_Series);
    tCPI_Series = tDataSample(:, [selForPloting{iCountry}, '_CPI']);
    tCPI_Series = rmmissing(tCPI_Series);
    
    nameCountry = tb_countries{find(strcmp(tb_countries.ISO2, selForPloting{iCountry})),'Name'}{:};
    if length(nameCountry) > 30
        nameAlt = tb_countries{find(strcmp(tb_countries.Name, nameCountry)), 'Name_Alternative_1'}{:};
        if length(nameAlt) <  length(nameCountry)
            nameCountry = nameAlt;
        end
    end
    
    tSampleDes(selForPloting{iCountry}, :) = ...
        { ...
        nameCountry, ...
        tb_countries{find(strcmp(tb_countries.ISO2, selForPloting{iCountry})),'Status'}, ...        
        datestr(tNR_Series.Time(1),'yyyyQQ'), ...
        datestr(tNR_Series.Time(end), 'yyyyQQ'), ...
        datestr(tCPI_Series.Time(1), 'yyyyQQ'), ...
        datestr(tCPI_Series.Time(end), 'yyyyQQ')};
end
tSampleDes
sSampleDes                  = struct();
sSampleDes.data             = cellfun( @(x) char(x), table2cell(tSampleDes) , 'UniformOutput', false);
sSampleDes.tableColLabels   = cellfun(@(x) strrep(x, '_', ' '), tSampleDes.Properties.VariableNames, 'UniformOutput', false);
sSampleDes.tableRowLabels   = tSampleDes.Properties.RowNames;
sSampleDes.tablePositioning = 'H';
sSampleDes.tableLabel       = 'tab_FiltersExercise_DataSample';
sSampleDes.tableCaption     = 'Filters: Summary of the data sample';
lTable                      = latexTable(sSampleDes);


%%%%%%%%%%%%%% Consolidate countries
countriesValidMdls = fieldnames(ts_All_Filters);
nValidMdls = length(countriesValidMdls);
ts_Consolidated = timetable();
for jj = 1:nValidMdls
    ts_aux = ts_All_Filters.(countriesValidMdls{jj})(:,'trend_Median');
    ts_aux.Properties.VariableNames = countriesValidMdls(jj);
    ts_Consolidated = synchronize(ts_Consolidated, ts_aux);
end

ts_Advanced = ts_Consolidated(:, selAdvancedCountries_2L( ...
    ismember(selAdvancedCountries_2L, ts_Consolidated.Properties.VariableNames)));
ts_Emerging = ts_Consolidated(:, selEmergingCountries_2L( ...
    ismember(selEmergingCountries_2L, ts_Consolidated.Properties.VariableNames)));
countriesValidMdls_Advanced = countriesValidMdls(ismember(countriesValidMdls, selAdvancedCountries_2L));
countriesValidMdls_Emerging = countriesValidMdls(ismember(countriesValidMdls, selEmergingCountries_2L));

r_median_Advanced           = median(ts_Advanced.Variables, 2, 'omitnan');
r_mean_Advanced             = mean(ts_Advanced.Variables, 2, 'omitnan');
r_weightedMean_Advanced     = pub_WeightedMean(ts_Advanced.Variables, ...
                    countryWeight_Advanced(countriesValidMdls_Advanced,:).Variables, 2);

r_median_Emerging           = median(ts_Emerging.Variables, 2, 'omitnan');
r_mean_Emerging             = mean(ts_Emerging.Variables, 2, 'omitnan');
r_weightedMean_Emerging     = pub_WeightedMean(ts_Emerging.Variables, ...
                    countryWeight_Emerging(countriesValidMdls_Emerging,:).Variables, 2);
                
ts_Consolidated = array2timetable([r_median_Advanced, r_mean_Advanced, r_weightedMean_Advanced, ...
    r_median_Emerging, r_mean_Emerging, r_weightedMean_Emerging], ...
    'RowTimes', ts_Consolidated.Properties.RowTimes, ...
    'VariableNames', {'advanced_median', 'advanced_mean', 'advanced_weightedMean', ...
    'emerging_median', 'emerging_mean', 'emerging_weightedMean'});


%%%%%%% Plot "Real neutral rate estimated by univariate filters"
f = figure;
previousInterpreter = pub_GraphSetInterpreter('latex');
p = plot(ts_Consolidated.Time, ts_Consolidated.Variables);
set(gca, 'FontSize', 16);
p(1).DisplayName = 'Advanced: Median of filters'' median';
p(2).DisplayName = 'Advanced: Mean of filters'' median';
p(3).DisplayName = 'Advanced: Mean of filters'' median weighted by GDP';
p(1).LineStyle = '-';
p(2).LineStyle = ':';
p(3).LineStyle = '--';
p(1).Color = 'blue';
p(2).Color = 'blue';
p(3).Color = 'blue';
p(4).DisplayName = 'Emerging: Median of filters'' median';
p(5).DisplayName = 'Emerging: Mean of filters'' median';
p(6).DisplayName = 'Emerging: Mean of filters'' median weighted by GDP';
p(4).LineStyle = '-';
p(5).LineStyle = ':';
p(6).LineStyle = '--';
p(4).Color = 'red';
p(5).Color = 'red';
p(6).Color = 'red';
p(4).LineWidth = 2;
p(5).LineWidth = 2;
p(6).LineWidth = 2;
xlim([datetime(2000,1,1) datetime(2020,1,1)]);
legend('location', 'northeast', 'Interpreter', 'latex');
pub_GraphDrawZeroAxis(p);
xlabel('date', 'Interpreter', 'latex');
ylabel('\% per annum', 'Interpreter', 'latex');
title('Real neutral rate estimated by univariate filters', 'Interpreter', 'latex');
graphName = 'rStar_Filters.png';
set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, graphName]);
pub_GraphSetInterpreter(previousInterpreter);


%% Test the difference between Advanced and Emerging

testSample = rmmissing(ts_Consolidated);
testSample = testSample(datetime(2000,1,1):datetime(2020,12,31),:);

medianDiff          =  testSample.emerging_median - testSample.advanced_median;
meanDiff            =  testSample.emerging_mean - testSample.advanced_mean;
weightedMeanDiff    =  testSample.emerging_weightedMean - testSample.advanced_weightedMean;

%%%%% Plot the differences (it is not in the paper)
f = figure;
p = plot(testSample.Time, [medianDiff, meanDiff, weightedMeanDiff]);
pub_GraphDrawZeroAxis(p);
legend('Median', 'Mean', 'Weighted Mean');
title('Difference between Emerging and Advanced real neutral rates');
ylabel('% per annum');
xlabel('date');

%%%% FIRST TEST (Comparison of means)

% THE 2 LINES BELOW MAY FIX A PATH PROBLEM WITH TTEST FUNCTION
% restoredefaultpath 
% rehash toolboxcache

% returns a test decision for the null hypothesis that the
% pairwise difference between data vectors
% x and y has a mean equal to zero.
% The result h is 1 if the test rejects the null hypothesis and 0 otherwise.
[h, p, ci, stats] = ttest( ...
                            testSample.emerging_median, ...
                            testSample.advanced_median, ...
                            'Tail', 'right', 'Alpha', 0.01);
disp('Test decision for the null hypothesis that the pairwise difference');
disp('between data vectors x and y has a mean equal to zero');                    
disp(h);     

% returns a test decision for the null hypothesis that the
% pairwise difference between data vectors
% x and y has a mean equal to zero.
% The result h is 1 if the test rejects the null hypothesis and 0 otherwise.
[h, p, ci, stats] = ttest( ...
                            testSample.emerging_mean, ...
                            testSample.advanced_mean, ...
                            'Tail', 'right', 'Alpha', 0.01);
disp('Test decision for the null hypothesis that the pairwise difference');
disp('between data vectors x and y has a mean equal to zero');                    
disp(h);

% returns a test decision for the null hypothesis that the
% pairwise difference between data vectors
% x and y has a mean equal to zero.
% The result h is 1 if the test rejects the null hypothesis and 0 otherwise.
[h, p, ci, stats] = ttest( ...
                            testSample.emerging_weightedMean, ...
                            testSample.advanced_weightedMean, ...
                            'Tail', 'right', 'Alpha', 0.01);
disp('Test decision for the null hypothesis that the pairwise difference');
disp('between data vectors x and y has a mean equal to zero');                    
disp(h);

%%%% SECOND TEST (ARIMA)

infoCriterion = 'bic';
maxP = 4;
maxD = 1;
maxQ = 4;

% Fit ARIMA on the difference (check significance of the constant)
disp('Median');
y           = medianDiff;
bestSpec    = pub_ARIMABestFit(y,maxP,maxD,maxQ,infoCriterion, 0,{},{});
mdl         = arima(bestSpec(1), bestSpec(2), bestSpec(3));
estMdl      = estimate(mdl, y);

% Fit ARIMA on the difference (check significance of the constant)
disp('Mean');
y           = meanDiff;
bestSpec    = pub_ARIMABestFit(y,maxP,maxD,maxQ,infoCriterion, 0,{},{});
mdl         = arima(bestSpec(1), bestSpec(2), bestSpec(3));
estMdl      = estimate(mdl, y);

% Fit ARIMA on the difference (check significance of the constant)
disp('Weighted Mean');
y           = weightedMeanDiff;
bestSpec    = pub_ARIMABestFit(y,maxP,maxD,maxQ,infoCriterion, 0,{},{});
mdl         = arima(bestSpec(1), bestSpec(2), bestSpec(3));
[estMdl,EstParamCov,logL,info]      = estimate(mdl, y);

%%%%%%%%%%%% SAVE
save([pathSaved, 'rStar_AllFilters.mat'], 'ts_All_Filters');