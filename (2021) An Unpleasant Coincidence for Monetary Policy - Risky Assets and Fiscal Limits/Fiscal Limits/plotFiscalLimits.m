%% Plot Fiscal Limits

plotActualHistograms = false;

if paperNumber == 2
    paper2_params_loadStandard;        % Parameter Values
elseif paperNumber == 3
    paper3_params_loadStandard;        % Parameter Values
end
params_loadValues;

approximateFiscalLimits;

%%% Set-up 1: steady-state         
ssValues_loadStandard;  % Steady-state values
ssValues_loadValues;
simulate_fiscalLimit;

% Plot histogram
f = figure;
set(gca,'LooseInset',get(gca,'TightInset'));
edu_GraphSetInterpreter('latex');
fiscalLimit_1 = fiscalLimit ./ (4 * ySS);
if plotActualHistograms
    h1 = histogram(fiscalLimit ./ (4 * ySS), 'Normalization', 'pdf', 'EdgeColor', 'none');
else
    pd = makedist('Normal', 'mu', fiscalLimit_mu./ (4 * ySS), 'sigma', fiscalLimit_std./ (4 * ySS));
    xx = [(fiscalLimit_mu - 5*fiscalLimit_std):0.001:(fiscalLimit_mu + 5*fiscalLimit_std)] ./ (4 * ySS);
    yy = pdf(pd, xx);
    h1 = area(xx*100, yy);
end
h1(1).LineStyle = '-';
h1(1).LineWidth = 1;
h1(1).FaceAlpha = 0.2;
%h1(1).Marker = 'none';
set(gca,'FontSize',14); % Scale fontsize of axes
title(['Fiscal limit distribution (', pub_ThousandSep(nSimulations, '%.0f', ','), ' simulations)'], 'interpreter', 'latex');
xlabel('max debt to steady-state annual output', 'interpreter', 'latex');
ylabel('probability density', 'interpreter', 'latex');
h1.DisplayName = [changeLabel, ' = ', num2str(simStruct.(changeType).(changeName))];
disp([changeLabel ' median of Plot 1:' num2str(median(fiscalLimit_1))]);
legend('Location', 'southoutside', 'Interpreter', 'latex', 'Orientation', 'horizontal');
grid('on');
%medLine = xline(median(h1.Data), '--k', num2str(round(median(h1.Data),2)), ...
%    'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left');
%medLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%%%%%%%%%%%

%%% Set-up 2: lower
simStruct.(changeType).(changeName) = paramLow;
eval(strcat(changeName, ' = ', string(paramLow), ';'));
ssValues_loadValues;
simulate_fiscalLimit;

% Plot histogram
hold on;
fiscalLimit_2 = fiscalLimit ./ (4 * ySS);
if plotActualHistograms
    h2 = histogram(fiscalLimit ./ (4 * ySS), 'Normalization', 'pdf', 'DisplayStyle', 'stairs');
else
    pd = makedist('Normal', 'mu', fiscalLimit_mu./ (4 * ySS), 'sigma', fiscalLimit_std./ (4 * ySS));
    xx = [(fiscalLimit_mu - 5*fiscalLimit_std):0.001:(fiscalLimit_mu + 5*fiscalLimit_std)] ./ (4 * ySS);
    yy = pdf(pd, xx);
    h2 = area(xx*100, yy);
end
h2(1).LineStyle = '--';
h2(1).LineWidth = 1;
h2(1).FaceAlpha = 0.2;
%h2(1).Marker = 'none';
h2.DisplayName = [changeLabel, ' = ', changeLabelLowValue];
disp([changeLabel ' median of Plot 2:' num2str(median(fiscalLimit_2))]);
%medLine = xline(median(h2.Data), '--k', num2str(round(median(h2.Data),2)), ...
%    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
%medLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%%%%%%%%%%%

%%% Set-up 3: greater
simStruct.(changeType).(changeName) = paramHigh;
eval(strcat(changeName, ' = ', string(paramHigh), ';'));
ssValues_loadValues;
simulate_fiscalLimit;

% Plot histogram
hold on;
fiscalLimit_3 = fiscalLimit ./ (4 * ySS);
if plotActualHistograms
    h3 = histogram(fiscalLimit ./ (4 * ySS), 'Normalization', 'pdf', 'DisplayStyle', 'stairs');
else
    pd = makedist('Normal', 'mu', fiscalLimit_mu./ (4 * ySS), 'sigma', fiscalLimit_std./ (4 * ySS));
    xx = [(fiscalLimit_mu - 5*fiscalLimit_std):0.001:(fiscalLimit_mu + 5*fiscalLimit_std)] ./ (4*ySS);
    yy = pdf(pd, xx);
    h3 = area(xx*100, yy);
end
h3(1).LineStyle = ':';
h3(1).LineWidth = 2;
h3(1).FaceAlpha = 0.2;
%h3(1).Marker = 'none';
h3.DisplayName = [changeLabel, ' = ',  changeLabelHighValue];
disp([changeLabel ' median of Plot 3:' num2str(median(fiscalLimit_3))]);
%medLine = xline(median(h3.Data), '--k', num2str(round(median(h3.Data),2)), ...
%    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
%medLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
%%%%%%%%%%%
set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

%%% Plot CDFs
f2 = figure;
set(gca,'LooseInset',get(gca,'TightInset'));
edu_GraphSetInterpreter('latex');
f21 = cdfplot(fiscalLimit_1);
f21 = plot(f21.XData*100, f21.YData*100);
set(gca,'FontSize',14); % Scale fontsize of axes
f21(1).LineStyle = '-';
f21(1).LineWidth = 1;
%f21(1).Marker = 'none';
hold on;
f22 = cdfplot(fiscalLimit_2);
f22 = plot(f22.XData*100, f22.YData*100);
f22(1).LineStyle = '--';
f22(1).LineWidth = 1;
%f22(1).Marker = 'none';
hold on;
f23 = cdfplot(fiscalLimit_3);
f23 = plot(f23.XData*100, f23.YData*100);
f23(1).LineStyle = ':';
f23(1).LineWidth = 2;
%f23(1).Marker = 'none';
xlim([50 250]);

title(['Fiscal limit cdf (', pub_ThousandSep(nSimulations, '%.0f', ','), ' simulations)'], 'interpreter', 'latex');
xlabel('max debt to steady-state annual output', 'interpreter', 'latex');
ylabel('cumm. prob. (\%)', 'interpreter', 'latex');

set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f2,[imagesFolder, simGraphName_CDF]);


%% Recover standard values

if paperNumber == 2
    paper2_params_loadStandard;        % Parameter Values
elseif paperNumber == 3
    paper3_params_loadStandard;        % Parameter Values
end
params_loadValues;

approximateFiscalLimits;
       
ssValues_loadStandard;  % Steady-state values
ssValues_loadValues;