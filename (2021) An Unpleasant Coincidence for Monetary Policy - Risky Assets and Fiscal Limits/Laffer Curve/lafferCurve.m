%% PAPER: An Unpleasant Coincidence for Monetary Policy: Risky Assets and Fiscal Limits
%
% Author: Eduardo G. C. Amaral
% Version date: 27/08/2020
%
% This routine plots the Laffer Curves

%% Plot the Laffer Curve

% Final or Temporary
saveToFinalFolder = false;

if saveToFinalFolder
    imagesFolder = pub_Path('/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 2/Images', 'C:\');
else
    imagesFolder = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) An Unpleasant Coincidence for Monetary Policy - Risky Assets and Fiscal Limits/Images', 'C:\');
end

f = figure;
pub_GraphSetInterpreter('latex');
tl = tiledlayout(1,6, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Set-up 1
eta     = 1;
kk      = 1;
aa      = 1;

expTax = @(tau, chi) tau * (1 - tau)^(1/chi) * (1/eta)^(1/chi) * (kk * aa)^(1+1/chi) ;
expTax_chi_050 = @(tau) expTax(tau, 0.5);
expTax_chi_075 = @(tau) expTax(tau, 0.75);
expTax_chi_100 = @(tau) expTax(tau, 1);

graph_1 = nexttile([1,2]);
fplot(expTax_chi_050, [0,1], 'LineWidth', 2);
hold on;
fplot(expTax_chi_075, [0,1], 'LineStyle', '--', 'LineWidth', 2);
hold on;
fplot(expTax_chi_100, [0,1], 'LineStyle', '-.', 'LineWidth', 2);

axLaffer = graph_1.Children;
for iLine = 1:length(axLaffer)
    [yMax, idxMax] = max(axLaffer(iLine).YData);
    xMax = axLaffer(iLine).XData(idxMax);
    linMax = xline(xMax, ':');
    linMax.Color = axLaffer(iLine).Color;
    linMax.LineWidth = 1.5;
end

set(graph_1, 'FontSize', 14);
title({'($\eta$ = 1.0; $A_t$ = 1.0)'}, 'FontSize', 14);
ylabel('$T$', 'FontSize', 14)
xlabel('$\tau$', 'FontSize', 14);

%% Set-up 2
eta     = 0.8;
kk      = 1;
aa  = 1;

expTax = @(tau, chi) tau * (1 - tau)^(1/chi) * (1/eta)^(1/chi) * (kk * aa)^(1+1/chi) ;
expTax_chi_050 = @(tau) expTax(tau, 0.5);
expTax_chi_075 = @(tau) expTax(tau, 0.75);
expTax_chi_100 = @(tau) expTax(tau, 1);

graph_2 = nexttile([1,2]);
fplot(expTax_chi_050, [0,1], 'LineWidth', 2);
hold on;
fplot(expTax_chi_075, [0,1], 'LineStyle', '--', 'LineWidth', 2);
hold on;
fplot(expTax_chi_100, [0,1], 'LineStyle', '-.', 'LineWidth', 2);

axLaffer = graph_2.Children;
for iLine = 1:length(axLaffer)
    [yMax, idxMax] = max(axLaffer(iLine).YData);
    xMax = axLaffer(iLine).XData(idxMax);
    linMax = xline(xMax, ':');
    linMax.Color = axLaffer(iLine).Color;
    linMax.LineWidth = 1.5;
end

set(graph_2, 'FontSize', 14);
title({'($\eta$ = 0.8; $A_t$ = 1.0)'}, 'FontSize', 14);
ylabel('$T$', 'FontSize', 14)
xlabel('$\tau$', 'FontSize', 14);

%% Set-up 3
eta     = 1;
kk      = 1;
aa  = 0.8;

expTax = @(tau, chi) tau * (1 - tau)^(1/chi) * (1/eta)^(1/chi) * (kk * aa)^(1+1/chi) ;
expTax_chi_050 = @(tau) expTax(tau, 0.5);
expTax_chi_075 = @(tau) expTax(tau, 0.75);
expTax_chi_100 = @(tau) expTax(tau, 1);

graph_3 = nexttile([1,2]);
fplot(expTax_chi_050, [0,1], 'LineWidth', 2);
hold on;
fplot(expTax_chi_075, [0,1], 'LineStyle', '--', 'LineWidth', 2);
hold on;
fplot(expTax_chi_100, [0,1], 'LineStyle', '-.', 'LineWidth', 2);

axLaffer = graph_3.Children;
for iLine = 1:length(axLaffer)
    [yMax, idxMax] = max(axLaffer(iLine).YData);
    xMax = axLaffer(iLine).XData(idxMax);
    linMax = xline(xMax, ':');
    linMax.Color = axLaffer(iLine).Color;
    linMax.LineWidth = 1.5;
end

set(graph_3, 'FontSize', 14);
title({'($\eta$ = 1.0; $A_t$ = 0.8)'}, 'FontSize', 14);
ylabel('$T$', 'FontSize', 14)
xlabel('$\tau$', 'FontSize', 14);

%% Final

% Link axes
linkaxes(tl.Children, 'xy');

% Add single legend slot
%poshL = get(graphLegend,'position');     % Getting its position
lgd = legend(...
             findobj(graph_1.Children, 'Type', 'FunctionLine'), ...
             {'$\chi$ = 1.00', '$\chi$ = 0.75', '$\chi$ = 0.50'}, ...
        'FontSize', 14, ...
        'Interpreter', 'latex', ...
        'Location', 'northeast', ...
        'Orientation', 'horizontal');
lgd.Layout.Tile = 'North';
lgd.Box = 'off';

set(f, 'Position',  [100, 100, 1000, 400]); % resize figure
exportgraphics(f, ...
    strjoin({imagesFolder 'Laffer curve' 'Laffer curves.png'}, filesep));