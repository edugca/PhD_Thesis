%% Load data

clear all; clc;

pathData    = 'Filters/Saved/';
pathImages  = pub_Path('Taylor Coefficient/Images/');

load([pathData, 'rStar_workspace.mat']);

countriesList = pub_List_Countries;

%% List countries

namesALL = ts.Properties.VariableNames;
tsCPI = ts(:, namesALL(startsWith(namesALL,'CPI_')));
tsNR = ts(:, namesALL(startsWith(namesALL,'NR_')));
tsY = ts(:, namesALL(startsWith(namesALL,'NY_')));

countriesCPI = cellfun(@(x) x(5:6), tsCPI.Properties.VariableNames, 'UniformOutput', false);
countriesNR = cellfun(@(x) x(4:5), tsNR.Properties.VariableNames, 'UniformOutput', false);

%% Find Taylor pii-coefficient

dateStart   = datetime(2000, 1, 1);
dateEnd     = datetime(2019, 10, 1);
sampleDates = dateStart:dateEnd;

responseVar = 'inflation'; % 'inflation', 'expectedInflation'

fontSize = 12;
fontSize_title = 12;
fontSize_labels = 10;
fontSize_legend = 12;

% Countries to remove from the sample (i.e. outliers)
remFromSample = {'TR'};

rulesList = {'noSmoothing', 'smoothing', 'trend', 'trendAndSmoothing'}; % 'noSmoothing', 'smoothing', 'trend', 'trendAndSmoothing'
rulesTitles = {'No trend nor smoothing', 'Smoothing', 'Trend', 'Trend and smoothing'};

if strcmp(responseVar, 'inflation')
    rulesSpecs = {
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \phi^{\pi}\left(\pi_t - \overline{\pi} \right)$',
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \alpha^ii_{t-1} + \left(1 - \alpha^i\right) \phi^{\pi}\left(\pi_t - \overline{\pi} \right)$',
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \alpha^tt + \phi^{\pi}\left(\pi_t - \overline{\pi} \right)$',
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \alpha^tt + \alpha^ii_{t-1} + \left(1 - \alpha^i\right) \phi^{\pi}\left(\pi_t - \overline{\pi} \right)$'
                   };
elseif strcmp(responseVar, 'expectedInflation')
    rulesSpecs = {
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \phi^{\pi}\left(\pi^e_t - \overline{\pi} \right)$',
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \alpha^ii_{t-1} + \left(1 - \alpha^i\right) \phi^{\pi}\left(\pi^e_t - \overline{\pi} \right)$',
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \alpha^tt + \phi^{\pi}\left(\pi^e_t - \overline{\pi} \right)$',
                    '$i_t = \left(\overline{r}+\overline{\pi}\right) + \alpha^tt + \alpha^ii_{t-1} + \left(1 - \alpha^i\right) \phi^{\pi}\left(\pi^e_t - \overline{\pi} \right)$'
                   };
end

f = figure;
nLins = 2;
nCols = 2;
tl = tiledlayout(nLins, nCols);
tl.TileSpacing  = 'compact';
tl.Padding      = 'compact';
selCountries    = {};
for iRule = 1:length(rulesList)

    taylorRule = rulesList{iRule};
    
    infCoeff   = NaN(length(countriesCPI),3);
    ctrStatus  = cell(length(countriesCPI),3); 
    for iCtr = 1:length(countriesCPI)
       ctrName = countriesCPI{iCtr};
       if ismember(ctrName, countriesNR) && ~ismember(ctrName, remFromSample)
           
           disp(['Running... Rule: ', taylorRule, '; Country: ', ctrName]);
           
           if iRule == 1
                selCountries = [selCountries, {ctrName}];
           end
           
            if strcmp(responseVar, 'inflation')
                responseVarName = ['CPI_' ctrName];
                responseVarNameVol = ['CPI_' ctrName];
            elseif strcmp(responseVar, 'expectedInflation')
                responseVarName = ['CPI_' ctrName '_Forecast_1Y'];
                responseVarNameVol = ['CPI_' ctrName];
            end
           
            %x = tsCPI{sampleDates,responseVarName} - movavg(tsCPI{sampleDates,responseVarName}, 'linear', 8);
            %y = tsNR{sampleDates,['NR_' ctrName]} - movavg(tsNR{sampleDates,['NR_' ctrName]}, 'linear', 8);
            %x = tsCPI{sampleDates,responseVarName} - mean(tsCPI{sampleDates,responseVarName}, 'omitnan');
            %y = tsNR{sampleDates,['NR_' ctrName]} - mean(tsNR{sampleDates,['NR_' ctrName]}, 'omitnan');
            x = (movprod((1 + tsCPI{sampleDates,responseVarName}./100).^(1/4), 4) - 1).*100;
            y = tsNR{sampleDates,['NR_' ctrName]};

            if strcmp(taylorRule, 'noSmoothing')
                mdl = fitlm(x - mean(x), y);
                infCoeff(iCtr,1) = mdl.Coefficients.Estimate(2);
                infCoeff(iCtr,2) = mdl.Coefficients.tStat(2);
            elseif strcmp(taylorRule, 'smoothing')
                mdl = fitlm([x(2:end) - mean(x), y(1:end-1)], y(2:end));
                infCoeff(iCtr,1) = mdl.Coefficients.Estimate(2) / (1 - mdl.Coefficients.Estimate(3));
                infCoeff(iCtr,2) = NaN;
            elseif strcmp(taylorRule, 'trend')
                mdl = fitlm([x - mean(x), (1:length(x))'], y);
                infCoeff(iCtr,1) = mdl.Coefficients.Estimate(2);
                infCoeff(iCtr,2) = mdl.Coefficients.tStat(2);
            elseif strcmp(taylorRule, 'trendAndSmoothing')
                mdl = fitlm([x(2:end) - mean(x), (2:length(x))', y(1:end-1)], y(2:end));
                infCoeff(iCtr,1) = mdl.Coefficients.Estimate(2) / (1 - mdl.Coefficients.Estimate(4));
                infCoeff(iCtr,2) = NaN;
            end

            infCoeff(iCtr,1) = mdl.Coefficients.Estimate(2);
            infCoeff(iCtr,2) = mdl.Coefficients.tStat(2);
            infCoeff(iCtr,3) = std(tsCPI{sampleDates,responseVarNameVol}, 'omitnan');

            ctrStatus{iCtr,1} = countriesList.ISO2{find(strcmp(ctrName, countriesList.ISO2))};
            ctrStatus{iCtr,2} = countriesList.Name{find(strcmp(ctrName, countriesList.ISO2))};
            ctrStatus{iCtr,3} = countriesList.Status{find(strcmp(ctrName, countriesList.ISO2))};
       end
    end

    %%%%%%% Plot results

    dx = 0.01;
    dy = 0.01;

    emeInfCoeff = infCoeff(strcmp(ctrStatus(:,3), 'Emerging'), :);
    advInfCoeff = infCoeff(strcmp(ctrStatus(:,3), 'Advanced'), :);

    emeCtrStatus = ctrStatus(strcmp(ctrStatus(:,3), 'Emerging'), :);
    advCtrStatus = ctrStatus(strcmp(ctrStatus(:,3), 'Advanced'), :);

    nexttile();
    pub_GraphSetInterpreter('latex');
    x = advInfCoeff(:,1);
    y = advInfCoeff(:,3);
    s1 = scatter(x, y, 'b', '*', 'DisplayName', 'Advanced economies');
    text(x+dx, y+dy, advCtrStatus(:,1), 'FontSize', fontSize_labels);
    hold on;
    x = emeInfCoeff(:,1);
    y = emeInfCoeff(:,3);
    s2 = scatter(x, y, 'r', 'x', 'DisplayName', 'Emerging economies');
    text(x+dx, y+dy, emeCtrStatus(:,1), 'FontSize', fontSize_labels);
    ls = lsline;
    ls(2).Color = [0 0 1];
    ls(1).Color = [1 0 0];
    ls(2).DisplayName = 'Advanced best-fit';
    ls(1).DisplayName = 'Emerging best-fit';
    uistack(ls(1), 'top');
    xlabel('$\phi^{\pi}$');
    ylabel('Inflation standard-deviation');
    set(gca, 'FontSize', fontSize);
    title([rulesTitles{iRule} newline rulesSpecs{iRule}], 'FontSize', fontSize_title);
    
    if iRule == 1
        lg = legend('FontSize', fontSize_legend, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'north';
        lg.Box = 'off';
    end
    
end

simGraphName = ['taylorRules_empirical_' responseVar '.png'];
set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
saveas(f,[pathImages, simGraphName]);