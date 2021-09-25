%%% Paper: Monetary Policy with Risky Assets
%%%
%%% Author: Eduardo G. C. Amaral
%%% Creation date: November 18, 2020
%%% Last modification: November 18, 2020
%%%
%%% Empirical evidence: Test the Power of Monetary Policy
%%%
%%% In case you find any mistake in this code, do not hesitate in sending
%%% me the heads up :)
%%%

clear all; clc;

imagesFolder = 'Empirical Correlations/Images/';

%% Load data

ts_rStar  = load('Filters/Saved/rStar_workspace.mat');

%data_LCRisk = load('/Users/Eduardo/OneDrive/MATLAB/My Thesis/Paper 5/Empirical Evidence/Paper Routines/determinantsLCRisk.mat');
%ts_riskLC = table2timetable(data_LCRisk.t(:,{'Date', 'defLC'}), 'RowTimes', 'Date');

ts_CDS = pub_Database_CDS_Fetch('Bloomberg');
ts_CDS.Properties.VariableNames = strcat('CDS_ ', ts_CDS.Properties.VariableNames);

ts_LCCS = pub_Database_LCCS_Fetch();
ts_LCCS.Properties.VariableNames = strcat('LCCS_ ', ts_LCCS.Properties.VariableNames);

load('Data/isoCountryCodes.mat', 't_ISOCountryCodes');

%% Save data

%load('dataset_empiricalCorrel.mat');
%save('dataset_empiricalCorrel.mat');

%% Organize data

% Country list
tb_countries = ts_rStar.tb_countries;
tb_countries.ShortName = tb_countries.Name;
tb_countries.Properties.RowNames = tb_countries.Name;
tb_countries{'United Kingdom of Great Britain and Northern Ireland', 'ShortName'} = {'United Kingdom'};
tb_countries{'United States of America', 'ShortName'} = {'United States'};

countryNames_CDS = cellfun(@(x) x(5:end), string(ts_CDS.Properties.VariableNames), 'UniformOutput', false);

advanced_conversion = ...
    {'US',  'UnitedStates';...
    'JP',   'Japan'; ...
    'DE',   'Germany'; ...
    'GB',   'UnitedKingdom'; ...
    'FR',   'France'; ...
    'IT',   'Italy'; ...
    'CA',   'Canada'; ...
    'KR',   'SouthKorea';...
    'ES',   'Spain'; ...
    'AU',   'Australia'; ...
    'NL',   'Netherlands'; ...
    'CH',   'Switzerland'; ...
    'TW',   'Taiwan'; ...
    'SE',   'Sweden'; ...
    'BE',   'Belgium'; ...
    'AT',   'Austria'; ...
    'NO',   'Norway'; ...
    'IL',   'Israel'; ...
    'IE',   'Ireland'; ...
    'HK',   'HongKong'};

emerging_conversion = ...
    {'CN', 'China'                    ; ...
    'IN', 'India'                     ; ...
    'BR', 'Brazil'                    ; ...
    'RU', 'Russia'                    ; ...
    'MX', 'Mexico'                    ; ...
    'ID', 'Indonesia'                 ; ...
    'SA', 'SaudiArabia'               ; ...
    'TR', 'Turkey'                    ; ...
    'PL', 'Poland'                    ; ...
    'TH', 'Thailand'                  ; ...
    'IR', 'Iran (Islamic Republic of)'; ...
    'NG', 'Nigeria'                   ; ...
    'AR', 'Argentina'                 ; ...
    'AE', 'United Arab Emirates'      ; ...
    'MY', 'Malaysia'                  ; ...
    'ZA', 'SouthAfrica'               ; ...
    'PH', 'Philippines'               ; ...
    'CO', 'Colombia'                  ; ...
    'BD', 'Bangladesh'                ; ...
    'EG', 'Egypt'                     };
    

country_conversion = [advanced_conversion; emerging_conversion];

selectedCountries = {};
for iCountry = 1:length(ts_CDS.Properties.VariableNames)
    cdsCountry = ts_CDS.Properties.VariableNames{iCountry}(5:end);
    if  matches(cdsCountry, country_conversion(:,2)) ~= 0
        countryCode = country_conversion{find(matches(country_conversion(:,2), cdsCountry)), 1};
        selectedCountries = [selectedCountries; {ts_CDS.Properties.VariableNames{iCountry}, countryCode}];
    end
end

for iCountry = 1:size(selectedCountries, 1)
    idx = find(matches(ts_CDS.Properties.VariableNames, selectedCountries(iCountry,1)));
    ts_CDS.Properties.VariableNames{idx} = ['CDS_', selectedCountries{iCountry,2}];
end

% CPI and Interest rate
cpiIdx = startsWith(ts_rStar.ts.Properties.VariableNames, 'CPI_').*(~endsWith(ts_rStar.ts.Properties.VariableNames, '_Forecast_1Y'));
nrIdx = startsWith(ts_rStar.ts.Properties.VariableNames, 'NR_');
tbIdx = startsWith(ts_rStar.ts.Properties.VariableNames, 'TB_');
fxavgIdx = startsWith(ts_rStar.ts.Properties.VariableNames, 'FXAVG_');
dateStart = datetime(2000, 1, 1);
dateEnd = datetime(2019, 10, 1);
ts_workData = ts_rStar.ts(dateStart:dateEnd,ts_rStar.ts.Properties.VariableNames(find(cpiIdx + nrIdx + tbIdx + fxavgIdx)));

% Risk premium
tbIdx = find(startsWith(ts_workData.Properties.VariableNames, 'TB_'));
usIdx = find(strcmp(ts_workData.Properties.VariableNames, 'TB_US'));
for ii = 1:length(tbIdx)
    ts_workData.(['SPREAD_', ts_workData.Properties.VariableNames{tbIdx(ii)}(4:5)]) = ts_workData{:, tbIdx(ii)} - ts_workData{:, usIdx};
end

% LC Default
%ts_workData = synchronize(ts_workData, retime(ts_riskLC, ts_workData.Properties.RowTimes));

% CDS (aggregate quarterly by mean)
ts_CDS_aux = retime(ts_CDS, ts_workData.Properties.RowTimes, 'mean');
ts_CDS_aux.Variables = ts_CDS_aux.Variables ./ 100;
ts_workData = synchronize(ts_workData, ts_CDS_aux);

% LCCS (aggregate quarterly by mean)
ts_LCCS_aux = retime(ts_LCCS, ts_workData.Properties.RowTimes, 'mean');
ts_LCCS_aux.Variables = ts_LCCS_aux.Variables ./ 100;
ts_workData = synchronize(ts_workData, ts_LCCS_aux);

% Probability of Default in LC
%ts_workData.defLC_BR = ts_workData.defLC * 100;
%ts_workData.defLC_BR_ABS1Q = ts_workData.defLC_BR - lagmatrix(ts_workData.defLC_BR,1);
%ts_workData.defLC_BR_4Q_ABS4Q = lagmatrix(ts_workData.defLC_BR,-4) - ts_workData.defLC_BR;

%% Scatter Countries

removeOutliers = true;
removeOutliersMethod = 'median';

% Configure plotting
%gap = [0.1 0.05];
%width_h = [0.1 0.05];
%width_w = [0.1 0.01];
%opt = {gap, width_h, width_w};
%subplot = @(m,n,p) subtightplot(m,n,p,opt{:});

allSeries = ts_workData.Properties.VariableNames;

countryGroupNames = {'Emerging', 'Advanced'};
countryGroupConversions = {emerging_conversion, advanced_conversion};

countryGroupMeasureInterestRate = {'NR', 'NR'};
countryGroupMeasureInflation    = {'CPI', 'CPI'};
countryGroupMeasureRisk         = {'LCCS', 'LCCS'}; %SPREAD, CDS or LCCS
countryGroupMeasureFX           = {'FXAVG', 'FXAVG'}; %FXAVG

countriesSelectedPerGroup = {};

tCorrelations = struct();
for iGroup = 1:length(countryGroupNames)
    
    groupName = countryGroupNames{iGroup};
    selected_conversion = countryGroupConversions{iGroup};
    countriesSelected = {};
    
    tCorrelations.(groupName) = cell2table(selected_conversion, 'VariableNames', {'ISO2', 'Name'});
    tCorrelations.(groupName).NR                = repmat("xxx", size(selected_conversion,1), 1);
    tCorrelations.(groupName).Pii               = repmat("xxx", size(selected_conversion,1), 1);
    tCorrelations.(groupName).Risk              = repmat("xxx", size(selected_conversion,1), 1);
    tCorrelations.(groupName).MinPoints_NR      = repmat("xxx", size(selected_conversion,1), 1);
    tCorrelations.(groupName).MinPoints_Pii     = repmat("xxx", size(selected_conversion,1), 1);
    tCorrelations.(groupName).MinPoints_Risk    = repmat("xxx", size(selected_conversion,1), 1);
    tCorrelations.(groupName).Corr_NR_Pii       = NaN(size(selected_conversion,1), 1);
    tCorrelations.(groupName).Corr_NR_Risk      = NaN(size(selected_conversion,1), 1);
    tCorrelations.(groupName).Corr_Pii_Risk     = NaN(size(selected_conversion,1), 1);
    tCorrelations.(groupName).Corr_Pii_Risk     = NaN(size(selected_conversion,1), 1);
    
    for iCountry = 1:size(selected_conversion,1)
        countryCode = selected_conversion{iCountry,1};

        codeNR = countryGroupMeasureInterestRate{iGroup};
        codePii = countryGroupMeasureInflation{iGroup};
        codeRisk = countryGroupMeasureRisk{iGroup};
        codeFX = countryGroupMeasureFX{iGroup};
        
        check_NR  = ismember([codeNR '_' countryCode], allSeries);
        check_Pii = ismember([codePii '_' countryCode], allSeries);
        check_Risk = ismember([codeRisk '_' countryCode], allSeries);
        
        % Minimum datapoints
        minObs = 10;
        if check_Risk == 1
            check_Risk  = sum(~isnan(ts_workData.([codeRisk '_' countryCode]))) > minObs;
        end
        if check_Pii == 1
            check_Pii  = sum(~isnan(ts_workData.([codePii '_' countryCode]))) > minObs;
        end
        if check_NR == 1
            check_NR  = sum(~isnan(ts_workData.([codeNR '_' countryCode]))) > minObs;
        end

        if check_NR && check_Pii && check_Risk
           countriesSelected = [countriesSelected; countryCode];
        end 
        
        % Fill in table
        countryIdx = find(strcmp(t_ISOCountryCodes{:,'Alpha2code'}, tCorrelations.(groupName){iCountry, 'ISO2'}));
        tCorrelations.(groupName){iCountry, 'Name'} = {t_ISOCountryCodes{countryIdx, 'Country'}};
        tCorrelations.(groupName){iCountry, 'NR'} = string(codeNR);
        tCorrelations.(groupName){iCountry, 'Pii'} = string(codePii);
        tCorrelations.(groupName){iCountry, 'Risk'} = string(codeRisk);
        tCorrelations.(groupName){iCountry, 'MinPoints_NR'} = check_NR;
        tCorrelations.(groupName){iCountry, 'MinPoints_Pii'} = check_Pii;
        tCorrelations.(groupName){iCountry, 'MinPoints_Risk'} = check_Risk;
        
    end

    countriesSelectedPerGroup{iGroup} = countriesSelected;
    
    % Plot
    iFigure = 0;
    nMaxLins = 4;
    fontSize = 10;
    previousInterpreter = pub_GraphSetInterpreter('Latex');
    for iCountry = 1:length(countriesSelected)

        iGraph = 0;

        % Create new figure
        if mod(iCountry, nMaxLins) == 1
            f = figure;
            nCols = 3;
            nLins = nMaxLins;
            iLine = 0;

            % [left right top bottom]  1 = 100%
            offsetSpecs = [0.1, 0.02 , 0.05 , 0.05];
            % [width height] 1 = 100%
            scaleSpecs = [0.8, 0.75];
            pos = iosr.figures.subfigrid(nLins, nCols, offsetSpecs, scaleSpecs);
        end

        iLine = iLine + 1;

        % Country codes
        cCode = countriesSelected{iCountry};
        cName = tb_countries{strcmp(tb_countries.ISO2, cCode), 'ShortName'}{:};
        cLine = find(strcmp(tCorrelations.(groupName).ISO2, cCode));
        
        % Interest rate
        ts_workData.([codeNR '_' cCode '_ABS1Q']) = ts_workData.([codeNR '_' cCode]) - lagmatrix(ts_workData.([codeNR '_' cCode]),1);
        ts_workData.([codeNR '_' cCode '_4Q_ABS4Q']) = lagmatrix(ts_workData.([codeNR '_' cCode]),-4) - ts_workData.([codeNR '_' cCode]);

        % CPI inflation
        ts_workData.([codePii '_' cCode '_ABS1Q']) = ts_workData.([codePii '_' cCode]) - lagmatrix(ts_workData.([codePii '_' cCode]),1);
        ts_workData.([codePii '_' cCode '_4Q_ABS4Q']) = lagmatrix(ts_workData.([codePii '_' cCode]),-4) - ts_workData.([codePii '_' cCode]);

        % Risk
        ts_workData.([codeRisk '_' cCode '_ABS1Q']) = ts_workData.([codeRisk '_' cCode]) - lagmatrix(ts_workData.([codeRisk '_' cCode]),1);
        ts_workData.([codeRisk '_' cCode '_4Q_ABS4Q']) = lagmatrix(ts_workData.([codeRisk '_' cCode]),-4) - ts_workData.([codeRisk '_' cCode]);
        
        % FX
        ts_workData.([codeFX '_' cCode '_ABS1Q']) = ts_workData.([codeFX '_' cCode]) - lagmatrix(ts_workData.([codeFX '_' cCode]),1);
        ts_workData.([codeFX '_' cCode '_4Q_ABS4Q']) = lagmatrix(ts_workData.([codeFX '_' cCode]),-4) - ts_workData.([codeFX '_' cCode]);
        
        %%%%%%% NR vs. Pii
        iGraph = iGraph + 1;
        subplot('position',pos(iLine,:,iGraph));
        %subplot(nLins, nCols, iGraph);
        x = ts_workData.([codeNR '_' cCode '_4Q_ABS4Q']);
        y = ts_workData.([codePii '_' cCode '_4Q_ABS4Q']);
        if removeOutliers
            [~, idxOutliers] = rmoutliers([x, y], removeOutliersMethod);
            x(idxOutliers)  = nan;
            y(idxOutliers)  = nan;
        end
        sc1 = scatter(x,y, 'filled');
        title([cName, newline, '4-quarter change'], 'FontWeight', 'normal');
        xlabel('Policy rate (p.p.)');
        ylabel('YoY CPI inflation (p.p.)');
        hline = lsline();
        mdl = fitlm(get(sc1,'xdata')', get(sc1,'ydata')');
        coeffs = mdl.Coefficients.Estimate;
        tStat  = mdl.Coefficients.tStat(2);
        equationStr = ['y = ' num2str(coeffs(1), '%.2f') ' + ' num2str(coeffs(2)), '(', num2str(tStat, '%.2f'), ')', 'x'];
        legend(hline, equationStr, 'Location', 'southoutside');
        legend('boxoff');

        dataCorr = rmmissing([get(sc1,'ydata')' get(sc1,'xdata')']);
        if ~isempty(dataCorr)
            corrCoeff = corr(dataCorr);
            tCorrelations.(groupName){cLine, 'Corr_NR_Pii'} =  corrCoeff(1,2);
        end
        
        %%%%%%% NR vs. Risk
        iGraph = iGraph + 1;
        subplot('position',pos(iLine,:,iGraph));
        %subplot(nLins, nCols, iGraph);
        x = ts_workData.([codeNR '_' cCode '_4Q_ABS4Q']);
        y = ts_workData.([codeRisk '_' cCode '_4Q_ABS4Q']);
        if removeOutliers
            [~, idxOutliers] = rmoutliers([x, y], removeOutliersMethod);
            x(idxOutliers)  = nan;
            y(idxOutliers)  = nan;
        end
        sc2 = scatter(x,y, 'filled');
        title([cName, newline, '4-quarter change'], 'FontWeight', 'normal');
        xlabel('Policy rate (p.p.)');
        ylabel([countryGroupMeasureRisk{iGroup} ' (p.p.)']);
        hline = lsline();
        mdl = fitlm(get(sc2,'xdata')', get(sc2,'ydata')');
        coeffs = mdl.Coefficients.Estimate;
        tStat  = mdl.Coefficients.tStat(2);
        equationStr = ['y = ' num2str(coeffs(1), '%.2f') ' + ' num2str(coeffs(2)), '(', num2str(tStat, '%.2f'), ')', 'x'];
        legend(hline, equationStr, 'Location', 'southoutside');
        legend('boxoff');
        
        dataCorr = rmmissing([get(sc2,'ydata')' get(sc2,'xdata')']);
        if ~isempty(dataCorr)
            corrCoeff = corr(dataCorr);
            tCorrelations.(groupName){cLine, 'Corr_NR_Risk'} =  corrCoeff(1,2);
        end
        
        %%%%%%% Pii vs. Risk
        iGraph = iGraph + 1;
        subplot('position',pos(iLine,:,iGraph));
        %subplot(nLins, nCols, iGraph);
        x = ts_workData.([codePii '_' cCode '_4Q_ABS4Q']);
        y = ts_workData.([codeRisk '_' cCode '_4Q_ABS4Q']);
        if removeOutliers
            [~, idxOutliers] = rmoutliers([x, y], removeOutliersMethod);
            x(idxOutliers)  = nan;
            y(idxOutliers)  = nan;
        end
        sc3 = scatter(x,y, 'filled');
        title([cName, newline, '4-quarter change'], 'FontWeight', 'normal');
        xlabel('YoY CPI inflation (p.p.)');
        ylabel([countryGroupMeasureRisk{iGroup} ' (p.p.)']);
        hline = lsline();
        mdl = fitlm(get(sc3,'xdata')', get(sc3,'ydata')');
        coeffs = mdl.Coefficients.Estimate;
        tStat  = mdl.Coefficients.tStat(2);
        equationStr = ['y = ' num2str(coeffs(1), '%.2f') ' + ' num2str(coeffs(2)), '(', num2str(tStat, '%.2f'), ')', 'x'];
        legend(hline, equationStr, 'Location', 'southoutside');
        legend('boxoff');
        
        dataCorr = rmmissing([get(sc3,'ydata')' get(sc3,'xdata')']);
        if ~isempty(dataCorr)
            corrCoeff = corr(dataCorr);
            tCorrelations.(groupName){cLine, 'Corr_Pii_Risk'} =  corrCoeff(1,2);
        end
        
        % Format all graphs
        if mod(iCountry, nMaxLins) == 0 || iCountry == length(countriesSelected)
            iFigure = iFigure + 1;
            set(f.Children, 'FontSize', fontSize);
            graphName = ['scatter_', groupName, '_NR_Inflation_Default_', num2str(iFigure), '.png'];
            set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
            saveas(f,[imagesFolder, graphName]);
        end

    end
end
pub_GraphSetInterpreter(previousInterpreter);

%% Scatter pooled economies

removeOutliers = true;
removeOutliersMethod = 'median';

for iGroup = 1:length(countryGroupNames)
    
    f = figure;
    previousInterpreter = pub_GraphSetInterpreter('latex');
    
    % Plot
    nCols = 1;
    nLins = 3;
    iGraph = 0;

    % [left right top bottom]  1 = 100%
    offsetSpecs = [0.1, 0.02 , 0.05 , 0.05];
    % [width height] 1 = 100%
    scaleSpecs = [0.85, 0.8];
    pos = iosr.figures.subfigrid(nLins, nCols, offsetSpecs, scaleSpecs);
    
    data_nr_cpi = [];
    data_nr_cds = [];
    data_cpi_cds = [];
    
    countriesSelected = countriesSelectedPerGroup{iGroup};
    codeNR      = countryGroupMeasureInterestRate{iGroup};
    codePii     = countryGroupMeasureInflation{iGroup};
    codeRisk    = countryGroupMeasureRisk{iGroup};
    codeFX      = countryGroupMeasureFX{iGroup};
    
    for iCountry = 1:length(countriesSelected)
        
        % Country codes
        cCode = countriesSelected{iCountry};
        cName = tb_countries{strcmp(tb_countries.ISO2, cCode),1}{:};

        nrData  = ts_workData.([codeNR '_' cCode '_4Q_ABS4Q']);
        cpiData = ts_workData.([codePii '_' cCode '_4Q_ABS4Q']);
        cdsData = ts_workData.([codeRisk '_' cCode '_4Q_ABS4Q']);
        fxData  = ts_workData.([codeFX '_' cCode '_4Q_ABS4Q']);

        if removeOutliers
            [~, idxOutliers] = rmoutliers([nrData, cpiData, cdsData], removeOutliersMethod);
            nrData(idxOutliers)  = nan;
            cpiData(idxOutliers)  = nan;
            cdsData(idxOutliers)  = nan;
        end
        
        isIT = pub_List_IsInflationTargeter(cName, ts_workData.Properties.RowTimes);

        data_nr_cpi     = [data_nr_cpi; ...
            rmmissing([repmat(iCountry,size(ts_workData,1),1), (1:size(ts_workData,1))', ...
            isIT, ...
            nrData, cpiData])];
        data_nr_cds     = [data_nr_cds; ...
            rmmissing([repmat(iCountry,size(ts_workData,1),1), (1:size(ts_workData,1))', ...
            isIT, ...
            nrData, cdsData])];
        data_cpi_cds    = [data_cpi_cds; ...
            rmmissing([repmat(iCountry,size(ts_workData,1),1), (1:size(ts_workData,1))', ...
            isIT, ...
            cpiData, cdsData, fxData])];

    end

    % Plot
    robustReg = 'off'; % 'on' or 'off' (outliers are supposed to have been removed above)

    iGraph = iGraph + 1;
    subplot('position',pos(iGraph,:,1));
    sc1 = scatter(data_nr_cpi(:,4), data_nr_cpi(:,5), 'filled');
    title([countryGroupNames{iGroup} ' economies: 4-quarter change'], 'FontWeight', 'normal');
    xlabel('Policy rate (p.p.)');
    ylabel('YoY CPI inflation (p.p.)');
    hline = lsline();
    mdl_nr_cpi = fitlm(get(sc1,'xdata')', get(sc1,'ydata')', ...
        'RobustOpts', robustReg, 'VarNames', {'nr','cpi'});
    mdl2_nr_cpi = fitlm([data_nr_cpi(:,1), get(sc1,'xdata')'], get(sc1,'ydata')', ...
        'CategoricalVars', 1, 'RobustOpts', robustReg, 'VarNames', {'country', 'nr','cpi'});
    mdl3_nr_cpi = fitlm([data_nr_cpi(:,[1,3]), get(sc1,'xdata')'], get(sc1,'ydata')', ...
        'CategoricalVars', 1, 'RobustOpts', robustReg, 'VarNames', {'country', 'inflationTarget', 'nr','cpi'});
    coeffs  = mdl_nr_cpi.Coefficients.Estimate;
    tStat   = mdl_nr_cpi.Coefficients.tStat;
    coeffs2 = mdl2_nr_cpi.Coefficients.Estimate;
    tStat2  = mdl2_nr_cpi.Coefficients.tStat;
    coeffs3 = mdl3_nr_cpi.Coefficients.Estimate;
    tStat3  = mdl3_nr_cpi.Coefficients.tStat;
    legStr = {  ['y = ' num2str(coeffs(1), '%.2f') ' + ' num2str(coeffs(2), '%.2f'), '(', num2str(tStat(end), '%.2f'), ')', 'x', newline, ...
                 'y = ' num2str(coeffs2(1), '%.2f') ' + country + ' num2str(coeffs2(end), '%.2f'), '(', num2str(tStat2(end), '%.2f'), ')', 'x', newline, ...
                 'y = ' num2str(coeffs3(1), '%.2f') ' + country + ' num2str(coeffs3(3), '%.2f') '(', num2str(tStat3(3), '%.2f'), ')' 'infTarget + ' num2str(coeffs3(end), '%.2f'), '(', num2str(tStat3(end), '%.2f'), ')', 'x'
                 ]};
    legend(hline, legStr, 'Location', 'southoutside');

    iGraph = iGraph + 1;
    subplot('position',pos(iGraph,:,1));
    sc2 = scatter(data_nr_cds(:,4), data_nr_cds(:,5), 'filled');
    title([countryGroupNames{iGroup} ' economies: 4-quarter change'], 'FontWeight', 'normal');
    xlabel('Policy rate (p.p.)');
    ylabel([codeRisk ' (p.p.)']);
    hline = lsline();
    mdl_nr_cds = fitlm(get(sc2,'xdata')', get(sc2,'ydata')', ...
        'RobustOpts', robustReg, 'VarNames', {'nr','cpi'});
    mdl2_nr_cds = fitlm([data_nr_cds(:,1), get(sc2,'xdata')'], get(sc2,'ydata')', ...
        'CategoricalVars', 1, 'RobustOpts', robustReg, 'VarNames', {'country', 'nr','cpi'});
    mdl3_nr_cds = fitlm([data_nr_cds(:,[1,3]), get(sc2,'xdata')'], get(sc2,'ydata')', ...
        'CategoricalVars', 1, 'RobustOpts', robustReg, 'VarNames', {'country', 'inflationTarget', 'nr','cpi'});
    coeffs  = mdl_nr_cds.Coefficients.Estimate;
    tStat   = mdl_nr_cds.Coefficients.tStat;
    coeffs2 = mdl2_nr_cds.Coefficients.Estimate;
    tStat2  = mdl2_nr_cds.Coefficients.tStat;
    coeffs3 = mdl3_nr_cds.Coefficients.Estimate;
    tStat3  = mdl3_nr_cds.Coefficients.tStat;
    legStr = {  ['y = ' num2str(coeffs(1), '%.2f') ' + ' num2str(coeffs(2), '%.2f'), '(', num2str(tStat(end), '%.2f'), ')', 'x', newline, ...
                 'y = ' num2str(coeffs2(1), '%.2f') ' + country + ' num2str(coeffs2(end), '%.2f'), '(', num2str(tStat2(end), '%.2f'), ')', 'x', newline, ...
                 'y = ' num2str(coeffs3(1), '%.2f') ' + country + ' num2str(coeffs3(3), '%.2f') '(', num2str(tStat3(3), '%.2f'), ')' 'infTarget + ' num2str(coeffs3(end), '%.2f'), '(', num2str(tStat3(end), '%.2f'), ')', 'x'
                 ]};
    legend(hline, legStr, 'Location', 'southoutside');

    iGraph = iGraph + 1;
    subplot('position',pos(iGraph,:,1));
    sc3 = scatter(data_cpi_cds(:,4), data_cpi_cds(:,5), 'filled');
    title([countryGroupNames{iGroup} ' economies: 4-quarter change'], 'FontWeight', 'normal');
    xlabel('YoY CPI inflation (p.p.)');
    ylabel([codeRisk ' (p.p.)']);
    hline = lsline();
    mdl_cpi_cds = fitlm(get(sc3,'xdata')', get(sc3,'ydata')', ...
        'RobustOpts', robustReg, 'VarNames', {'nr','cpi'});
    mdl2_cpi_cds = fitlm([data_cpi_cds(:,1), get(sc3,'xdata')'], get(sc3,'ydata')', ...
        'CategoricalVars', 1, 'RobustOpts', robustReg, 'VarNames', {'country', 'nr','cpi'});
    mdl3_cpi_cds = fitlm([data_cpi_cds(:,[1,3]), get(sc3,'xdata')'], get(sc3,'ydata')', ...
        'CategoricalVars', 1, 'RobustOpts', robustReg, 'VarNames', {'country', 'inflationTarget', 'nr','cpi'});
    coeffs  = mdl_cpi_cds.Coefficients.Estimate;
    tStat   = mdl_cpi_cds.Coefficients.tStat;
    coeffs2 = mdl2_cpi_cds.Coefficients.Estimate;
    tStat2  = mdl2_cpi_cds.Coefficients.tStat;
    coeffs3 = mdl3_cpi_cds.Coefficients.Estimate;
    tStat3  = mdl3_cpi_cds.Coefficients.tStat;
    legStr = {  ['y = ' num2str(coeffs(1), '%.2f') ' + ' num2str(coeffs(2), '%.2f'), '(', num2str(tStat(end), '%.2f'), ')', 'x', newline, ...
                 'y = ' num2str(coeffs2(1), '%.2f') ' + country + ' num2str(coeffs2(end), '%.2f'), '(', num2str(tStat2(end), '%.2f'), ')', 'x', newline, ...
                 'y = ' num2str(coeffs3(1), '%.2f') ' + country + ' num2str(coeffs3(3), '%.2f') '(', num2str(tStat3(3), '%.2f'), ')' 'infTarget + ' num2str(coeffs3(end), '%.2f'), '(', num2str(tStat3(end), '%.2f'), ')', 'x'
                 ]};
    legend(hline, legStr, 'Location', 'southoutside');

    % Format all graphs
    if mod(iCountry, nMaxLins) == 0 || iCountry == length(countriesSelected)
        iFigure = iFigure + 1;
        set(f.Children, 'FontSize', fontSize);
        graphName = ['scatter_' countryGroupNames{iGroup} 'All_NR_Inflation_Default.png'];
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        saveas(f,[imagesFolder, graphName]);
    end
    pub_GraphSetInterpreter(previousInterpreter);

    % Examine robust regression
    if strcmp(robustReg, 'on')
        % Find the index of the outlier.
        % Examine the weight of the outlier in the robust fit.
        % Compare the outlier weight with the median weight.
        disp('Outlier weight vs. Median Weight');
        [~,outlier] = max(mdl_nr_cpi.Residuals.Raw);
        mdl_nr_cpi.Robust.Weights(outlier)
        median(mdl_nr_cpi.Robust.Weights)

        [~,outlier] = max(mdl_nr_cds.Residuals.Raw);
        mdl_nr_cds.Robust.Weights(outlier)
        median(mdl_nr_cds.Robust.Weights)

        [~,outlier] = max(mdl_cpi_cds.Residuals.Raw);
        mdl_cpi_cds.Robust.Weights(outlier)
        median(mdl_cpi_cds.Robust.Weights)
    end

%     %%%%%%%%%%
%     % The following tests require the installation of the MatLab add-on "Panel Data Toolbox for MATLAB"
%     % https://www.mathworks.com/matlabcentral/fileexchange/51320-panel-data-toolbox-for-matlab
%     % Álvarez, Inmaculada C.; Barbero, Javier and Zofío, José L, (2017) A Panel Data Toolbox for MATLAB. Journal of Statistical Software. Volume 76, Issue 6, pp 1-27. http://dx.doi.org/10.18637/jss.v076.i06

%     
%     % Organize panel data
%     id   = data_nr_cpi(:,1);
%     time = data_nr_cpi(:,2);
%     X    = data_nr_cpi(:,4);
%     y    = data_nr_cpi(:,5);
% 
%     % Breusch-Pagan test
%     % It tests whether the variance of the errors from a regression is
%     % dependent on the values of the independent variables. In that case,
%     % heteroskedasticity is present.
% 
%     estPooledOLS = ols(y, X);
%     test = bphettest(estPooledOLS);
%     testdisp(test);
% 
%     estPooledOLS_timeDummy = ols(y, id);
%     test = bphettest(estPooledOLS_timeDummy);
%     testdisp(test);
% 
%     % Fixed effects
%     method = 'fe';
%     est_fe = panel( id, time, y, X, method);
% 
%     % Random effects
%     method = 're';
%     est_re = panel( id, time, y, X, method);
% 
%     % Hausman test
%     % The Hausman test can be also used to differentiate between
%     % fixed effects model and random effects model in panel data.
%     % In this case, Random effects (RE) is preferred under the null hypothesis
%     % due to higher efficiency, while under the alternative Fixed effects (FE)
%     % is at least as consistent and thus preferred.
%     test = hausmantest(est_fe, est_re);
%     testdisp(test);

end
