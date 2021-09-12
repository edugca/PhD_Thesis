%% COMPARISON WITH BLANCHARD (2004)

%% Housekeeping
clear;
clc;
close all;

pathMain    = pub_Path('/Users/Eduardo/OneDrive/MATLAB/Resources/Papers/(2021) An Unpleasant Coincidence for Monetary Policy - Risky Assets and Fiscal Limits', 'C:\');
pathTables  = [pathMain filesep 'Tables'];
pathImages  = [pathMain filesep 'Images'];

%% Plot graphs

fontSize = 12;
legendSize = 14;
pColors = linspecer(3); 
customizeLines = struct();
customizeLines.custom_LineWidth = {1, 2, 1, 1, 1, 1, 1};
customizeLines.custom_LineStyle = {'-', ':', '-.', '--', '-', ':', '-.'};
customizeLines.custom_LineColor = {'blue', 'red', 'black', 'blue', 'blue', 'blue', 'blue'};
customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};

f = figure;
set(gca, 'FontSize', 14);
previousInterpreter = edu_GraphSetInterpreter('latex');

% Set layout and reduce empty spaces in the graph
nLins = 1;
nCols = 2;
tl = tiledlayout(nLins, nCols);
tl.TileSpacing  = 'none';
tl.Padding      = 'none';

% Parameters
eeta = 0.9; % between 0 and 1
rStar = 0.02;
llambda = 0.1;
tthetaStar = 1.5;
primSurplus = 0.01;
mmu = 0.9;

syms r eepsilon p debt

capitalFlowEq(r, eepsilon, p) = (1+r) - eepsilon^(eeta - 1) * (1+rStar) - (1-llambda)*tthetaStar*p;
capitalFlowEq_risky(r, eepsilon, p) = (1-p)*(1+r) - eepsilon^(eeta - 1) * (1+rStar) - (1-llambda)*tthetaStar*p;
probDefault = ( (1+r)/(1-p) + llambda*tthetaStar*p/(1-p) ) * ( mmu*eepsilon + (1-mmu))*debt - primSurplus;
probDefault_risky(r, eepsilon, debt, p) = ( (1+r) + llambda*tthetaStar*p/(1-p) ) * ( mmu*eepsilon + (1-mmu))*debt - primSurplus;

solve(probDefault == p, eepsilon)

pub_GraphSetInterpreter('latex');

% Blanchard's Figure 3
f = figure;
p = fplot(@(eepsilon) probDefault(0.01, eepsilon, 0.5, 0.1)^2,[0,1], 'b');
hold on;
p = fplot(@(eepsilon) probDefault(0.10, eepsilon, 0.5, 0.1)^2,[0,1], 'b:', 'LineWidth', 2);
hold on;
p = fplot(@(eepsilon) probDefault(0.01, eepsilon, 1.0, 0.1)^2,[0,1], 'r');
hold on;
p = fplot(@(eepsilon) probDefault(0.10, eepsilon, 1.0, 0.1)^2,[0,1], 'r:', 'LineWidth', 2);
ylabel('p');
xlabel('$\epsilon$');


f = figure;
p = fplot(@(eepsilon) probDefault_risky(0.01, eepsilon, 0.5, 0.1)^2,[0,1], 'b');
hold on;
p = fplot(@(eepsilon) probDefault_risky(0.10, eepsilon, 0.5, 0.1)^2,[0,1], 'b:', 'LineWidth', 2);
hold on;
p = fplot(@(eepsilon) probDefault_risky(0.01, eepsilon, 1.0, 0.1)^2,[0,1], 'r');
hold on;
p = fplot(@(eepsilon) probDefault_risky(0.10, eepsilon, 1.0, 0.1)^2,[0,1], 'r:', 'LineWidth', 2);
ylabel('p');
xlabel('$\epsilon$');


% Blanchard's Figure 4a
f = figure;
p = fplot(@(p) finverse(probDefault(0.01, eepsilon, 0.5, p)^2),[0,1], 'b');
hold on;
p = fplot(@(p) probDefault(0.10, eepsilon, 0.5, p)^2,[0,1], 'b:', 'LineWidth', 2);
hold on;
p = fplot(@(p) probDefault(0.01, eepsilon, 1.0, p)^2,[0,1], 'r');
hold on;
p = fplot(@(p) probDefault(0.10, eepsilon, 1.0, p)^2,[0,1], 'r:', 'LineWidth', 2);
ylabel('$\epsilon$');
xlabel('p');

%%
for iGraph = 1

    nexttile();
    
    % Plot the functions
    xMin = 0.4;
    xMax = 1.3;
    yMin = 0.5;
    yMax = 2.7;
    x = xMin:0.01:xMax;
    
    if iGraph == 1
        plot(x, fisherEq(1, x), 'k', ...
             x, fisherEq(SR_val, x), 'k--', ...
             x, riskAdjPolRule(1, x, withELB), 'r', ...
             'linewidth', 2);
    elseif iGraph == 2
        plot(x, fisherEq(1, x), 'k', ...
            x, fisherEq(SR_val, x), 'k--', ...
            x, riskAdjPolRule(1, x, withELB), 'r', ...
            x, riskAdjPolRule(SR_val, x, withELB), 'r--', ...
            'linewidth', 2);
    end
         
    % get rid of standard axes decorations
    set(gca, 'Xtick', [], 'Ytick', [], 'box', 'off');
    set(gca,'visible','off')
    
    % Fix the axes sizes
    axis([xMin-0.2 xMax+0.2 yMin-0.2 yMax+0.2]);

    % Point where they meet
    xDSS = fzero(@(x) fisherEq(1, x) - riskAdjPolRule(1, x, withELB),1);
    yDSS = fisherEq(1, xDSS);

    xEM = fzero(@(x) fisherEq(SR_val, x) - riskAdjPolRule(1, x, withELB),1);
    yEM = fisherEq(SR_val, xEM);
    
    xEM2 = fzero(@(x) fisherEq(SR_val, x) - riskAdjPolRule(SR_val, x, withELB),1);
    yEM2 = fisherEq(SR_val, xEM2);

    % the equilibrium point
    hold on;
    plot(xDSS, yDSS, 'k.', 'markersize', 20);
    hold on;
    plot(xEM, yEM, 'k.', 'markersize', 20)
    
    if iGraph ==2
        hold on;
        plot(xEM2, yEM2, 'k.', 'markersize', 20)
    end

    % the locational lines
    hold on;
    line([xDSS xMin; xDSS xDSS], [yMin yDSS; yDSS yDSS], 'linestyle', ':', 'color', 'k', 'LineWidth', 1.2);
    if iGraph == 1
        line([xEM xMin; xEM xEM], [yMin yEM; yEM yEM], 'linestyle', ':', 'color', 'k', 'LineWidth', 1.2);
    else
        line([xEM2 xMin; xEM2 xEM2], [yMin yEM2; yEM2 yEM2], 'linestyle', ':', 'color', 'k', 'LineWidth', 1.2);
    end

    % the arrows
    hold on;
    q = quiver(xMin, yMin, 0, yMax-yMin, 0, '-k', 'LineWidth', 2, 'MaxHeadSize', 0.2/(yMax-yMin), 'Clipping', 'off');
    hold on;
    q = quiver(xMin, yMin, xMax-xMin, 0, 0, '-k', 'LineWidth', 2, 'MaxHeadSize', 0.2/(yMax-yMin), 'Clipping', 'off');

    % the static texts
    if withELB
        text(xMin - 0.1, 1,            'ELB', 'fontweight', 'bold', 'FontSize', 12);
    end
    if iGraph == 1
        text(PiiTarg-0.008, yMin-0.1, '$\overline{\Pi}$', 'fontweight', 'bold', 'FontSize', 12);
        text(xEM-0.008, yMin-0.1, '$\Pi^{EM}$', 'fontweight', 'bold', 'FontSize', 12);
        text(xMin-0.1, yEM, '$R^{EM}$', 'fontweight', 'bold', 'FontSize', 12);
        text(xMin-0.05, yDSS, '$\overline{R}$', 'fontweight', 'bold', 'FontSize', 12);
    elseif iGraph == 2
        text(PiiTarg-0.008, yMin-0.1, '$\overline{\Pi} = \Pi^{EM''}$', 'fontweight', 'bold', 'FontSize', 12);
        text(xMin-0.1, yEM2, '$R^{EM''}$', 'fontweight', 'bold', 'FontSize', 12);
        text(xMin-0.05, yDSS, '$\overline{R}$', 'fontweight', 'bold', 'FontSize', 12);
    end
    text(xMax + 0.015, yMin,        '$\Pi$', 'fontweight', 'bold', 'FontSize', 14);
    text(xMin-0.02, yMax + 0.07,    '$R$', 'fontweight', 'bold', 'FontSize', 14);

    text(xDSS+0.02, yDSS-0.03,      'DSS', 'fontweight', 'bold', 'FontSize', 12);
    text(xEM+0.03, yEM-0.03,        'EM', 'fontweight', 'bold', 'FontSize', 12);
    
    if iGraph == 2
        text(xEM2-0.08, yEM2+0.06,       'EM''', 'fontweight', 'bold', 'FontSize', 12);
    end

    title('Fisher Relation and Taylor Rule', 'FontSize', 14);

end

lg = legend('Fisher Relation', 'Risk-Adjusted Fisher Relation', 'Taylor Rule', 'Risk-Adjusted Taylor Rule', ...
        'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', legendSize);
lg.Layout.Tile = 'North';
lg.Box = 'off';

% Export
set(f, 'Position',  [100, 0, 1000, 400]); % resize figure
exportgraphics(f, ...
    strjoin({pathImages 'Comparison' 'Graph - Stylized Taylor Rule.png'}, filesep));
edu_GraphSetInterpreter(previousInterpreter);