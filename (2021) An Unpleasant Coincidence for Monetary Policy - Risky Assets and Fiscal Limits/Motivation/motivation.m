%% Plot historical data

% Paths of the model
pathMain    = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) An Unpleasant Coincidence for Monetary Policy - Risky Assets and Fiscal Limits', 'C:\');
pathImages  = [pathMain filesep 'Images' filesep 'Motivation'];

% Read data file
pathHistorical = [pathMain filesep 'Motivation' filesep 'Historical data.xlsx'];
tHist = readtable(pathHistorical, 'HeaderLines', 2);

verticalLines = {
                    datetime(1986,3,1), 'Cruzado Plan';
                    datetime(1990,3,1), 'Collor Plan';
                    datetime(1994,6,1), 'Real Plan';
                  };

shadedAreas = {
                    datetime(1980,1,1), datetime(1986,3,1), [0.6510    0.6510    0.6510], 'Chronic inflation';
                    datetime(1986,3,1), datetime(1994,7,1), [0.8510    0.3255    0.0980], 'Economic plans';
                    datetime(1994,7,1), datetime(1999,6,1), [0    0.4471    0.7412], 'Fixed exchange rate';
                    datetime(1999,6,1), datetime(2021,6,1), [0.4667    0.6745    0.1882], 'Inflation targeting';
                  };

% Plot historical data
f = figure;
tl = tiledlayout(3, 1, 'TileSpacing', 'none', 'Padding', 'compact');
pub_GraphSetInterpreter('latex');
h = gobjects(4,1); 

%%%%%%%% GRAPH 1
nexttile();
yMax = 90;
yMin = -10;
p = plot(tHist.DATE, [tHist.POLICY_MOM tHist.CPI_MOM]);
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(2).LineStyle = '-';
set(gca,'FontSize',14);
title({'Nominal Interest Rate and Inflation'}, 'interpreter', 'latex');
ylabel('\% MoM', 'interpreter', 'latex');
ylim([yMin yMax]);
yticks(yMin:20:yMax);
xlim([datetime(1980,1,1) datetime(2020,1,1)]);
xtickformat('MMM-yy')
%xticks(datetime(1980,3,1):calmonths(36):datetime(2021,3,1));
ytickformat('%.0f');
%xtickformat('%.0f');
p(1).Parent.XGrid = 'on';
p(1).Parent.YGrid = 'on';
legend({'Overnight nominal interest rate', 'CPI inflation rate'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Color', 'none', 'FontSize', 14);
legend boxoff;
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
for iDate = 1:size(verticalLines, 1)
    xl = xline(verticalLines{iDate, 1}, ':k', verticalLines{iDate, 2}, 'HandleVisibility','off', 'FontSize', 12);
    xl.LineWidth = 1.5;
    xl.LabelVerticalAlignment = 'top';
    xl.LabelHorizontalAlignment = 'right';
    xl.Interpreter = 'latex';
    %xl.FontWeight = 'bold';
end
for iDate = 1:size(shadedAreas, 1)
    xl = fill(  [shadedAreas{iDate, 1} shadedAreas{iDate, 1} shadedAreas{iDate, 2} shadedAreas{iDate, 2}], ...
                [yMin yMax yMax yMin], '', 'LineStyle', 'none', 'FaceColor', shadedAreas{iDate, 3}, 'FaceAlpha', 0.15, ...
                'HandleVisibility', 'off', 'DisplayName', shadedAreas{iDate, 4});
end
hold off;
pub_GraphDrawZeroAxis(p(1));

%%%%%%%% GRAPH 2
nexttile();
yMax = 100;
yMin = -100;
p = plot(tHist.DATE, [tHist.REAL_QOQAN]);
p(1).LineWidth = 2;
p(1).Color = [0/255 102/255 51/255];
set(gca,'FontSize',14);
title({'Ex-Post Real Interest Rate'}, 'interpreter', 'latex');
ylabel('\% QoQ (annualized)', 'interpreter', 'latex');
ylim([yMin yMax]);
yticks(yMin:25:yMax);
xlim([datetime(1980,1,1) datetime(2020,1,1)]);
xtickformat('MMM-yy')
%xticks(datetime(1980,3,1):calmonths(36):datetime(2021,3,1));
ytickformat('%.0f');
%xtickformat('%.0f');
p(1).Parent.XGrid = 'on';
p(1).Parent.YGrid = 'on';
%legend({'Ex-post real interest rate'}, ...
%    'Location', 'northoutside', 'Orientation', 'horizontal', 'Color', 'none', 'FontSize', 14);
%legend boxoff;
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
for iDate = 1:size(verticalLines, 1)
    xl = xline(verticalLines{iDate, 1}, ':k', verticalLines{iDate, 2}, 'HandleVisibility','off', 'FontSize', 12);
    xl.LineWidth = 1.5;
    xl.LabelVerticalAlignment = 'top';
    xl.LabelHorizontalAlignment = 'right';
    xl.Interpreter = 'latex';
    %xl.FontWeight = 'bold';
end
for iDate = 1:size(shadedAreas, 1)
    xl = fill(  [shadedAreas{iDate, 1} shadedAreas{iDate, 1} shadedAreas{iDate, 2} shadedAreas{iDate, 2}], ...
                [yMin yMax yMax yMin], '', 'LineStyle', 'none', 'FaceColor', shadedAreas{iDate, 3}, 'FaceAlpha', 0.15, ...
                'HandleVisibility', 'off', 'DisplayName', shadedAreas{iDate, 4});        
    h(iDate) = xl;
end
hold off;
pub_GraphDrawZeroAxis(p(1));


%%%%%%%% GRAPH 3
nexttile();
yMax = 90;
yMin = -10;
p = plot(tHist.DATE, [tHist.POLICY_QOQAN, tHist.CPI_QOQAN]);
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(2).LineStyle = '-';
set(gca,'FontSize',14);
title({'Nominal Interest Rate and Inflation (after the Real Plan)'}, 'interpreter', 'latex');
ylabel('\% QoQ (annualized)', 'interpreter', 'latex');
ylim([yMin yMax]);
yticks(yMin:20:yMax);
xlim([datetime(1994,6,1) datetime(2020,1,1)]);
xtickformat('MMM-yy')
%xticks(datetime(1980,3,1):calmonths(36):datetime(2021,3,1));
ytickformat('%.0f');
%xtickformat('%.0f');
p(1).Parent.XGrid = 'on';
p(1).Parent.YGrid = 'on';
legend({'Overnight nominal interest rate', 'CPI inflation rate'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Color', 'none', 'FontSize', 14);
legend boxoff;
set(gca,'LooseInset',get(gca,'TightInset'));
hold on;
% for iDate = 1:size(verticalLines, 1)
%     xl = xline(verticalLines{iDate, 1}, ':k', verticalLines{iDate, 2}, 'HandleVisibility','off', 'FontSize', 12);
%     xl.LineWidth = 1.5;
%     xl.LabelVerticalAlignment = 'top';
%     xl.LabelHorizontalAlignment = 'right';
%     xl.Interpreter = 'latex';
%     %xl.FontWeight = 'bold';
% end
for iDate = 1:size(shadedAreas, 1)
    xl = fill(  [shadedAreas{iDate, 1} shadedAreas{iDate, 1} shadedAreas{iDate, 2} shadedAreas{iDate, 2}], ...
                [yMin yMax yMax yMin], '', 'LineStyle', 'none', 'FaceColor', shadedAreas{iDate, 3}, 'FaceAlpha', 0.15, ...
                'HandleVisibility', 'off', 'DisplayName', shadedAreas{iDate, 4});
end
hold off;
pub_GraphDrawZeroAxis(p(1));

lg  = legend(h(:), 'Orientation', 'horizontal', 'Color', 'none', 'FontSize', 14); 
lg.Layout.Tile = 'North'; % <-- place legend east of tiles
lg.Box = 'off';

%%%%%%%%%%%
simGraphName = 'historicalData_Brazil.png';
set(gcf, 'Position',  [100, 0, 1000, 900]); % resize figure
saveas(f,[pathImages filesep simGraphName]);
