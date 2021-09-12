%% Wicksellian model augmented with default in the natural rate

isFinal = false;

if isFinal
    imagesFolder = '/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Images/Partial Equilibrium/';
else
    imagesFolder = 'Model/Images/';
end

%%% To replicate the Paper's figures, change phi (or phiPii) and set the respective policyTarget rule
phi         = 0.1; % price-level targeting: monetary policy parameter > 0
phiPii      = 1.1; % inflation targeting:monetary policy parameter > 1

policyTargetRule = 'inflation'; % priceLevel, inflation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha       = 1/(1+phi);
alphaPii    = 1/phiPii;

nPeriods        = 100;
nFuture         = 100;
nSimulations    = 4;
p               = nan(nPeriods, nSimulations);
pii             = nan(nPeriods, nSimulations);
rN              = nan(nPeriods + nFuture,1);
iota            = nan(nPeriods + nFuture,1);
nrGap           = nan(nPeriods + nFuture,1);
defPol          = nan(nPeriods + nFuture,1);
nrBad           = nan(nPeriods + nFuture,1);
defR            = nan(nPeriods + nFuture,1);
rGood           = nan(nPeriods + nFuture,1);
rBad            = nan(nPeriods + nFuture,1);

fontSize        = 14;
titleFontSize   = 16;

%% Partial equilibrium with risky natural

%%% Steady state
iSimul       = 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end

   if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
   
end
%%%%%%%%

%%% Jump in defR
iSimul        = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            defR(j)        = 0.15;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in defR (MIT shocks)
iSimul        = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            defR(ii)        = 0.15;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end
    
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rGood
iSimul        = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            rGood(j)        = 0.06;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end
    
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rGood (MIT shocks)
iSimul        = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            rGood(ii)        = 0.06;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end

    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rBad
iSimul        = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            rBad(j)        = -0.50;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in rBad (MIT shocks)
iSimul        = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            rBad(ii)        = -0.50;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end
   
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in iota
iSimul       = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            iota(j)        = 0.05;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end
    
    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% Jump in iota (MIT shocks)
iSimul       = iSimul + 1;
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);
iota(:)      = rN(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            iota(ii)        = 0.05;
        end
        
        rN(j) = (1-defR(j))*rGood(j) + defR(j)*rBad(j);
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end

    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% rN stochastic; Mon pol does not track the natural interest rate
iSimul       = iSimul + 1;
iota(:)      = 0.04;
defR(:)      = max(0, 0.1 + 0.1*randn(length(rN),1));
rGood(:)     = 0.04;
rBad(:)      = -0.10;
rN(:)        = (1-defR(:)).*rGood(:) + defR(:).*rBad(:);

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end

    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end
%%%%%%%%

%%% rN stochastic; Mon pol tracks the natural interest rate
iSimul       = iSimul + 1;
rN(:)        = rN(:); % use the same sequence of shocks
iota(:)      = rN(:);
defR(:)      = 0.10;
rGood(:)     = 0.04;
rBad(:)      = -0.60;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(j) - iota(j));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(j) - iota(j));
        end
        
    end

    if strcmpi(policyTargetRule, 'priceLevel')
        p(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            pii(ii, iSimul) =  p(ii, iSimul) - 0;
        else
            pii(ii, iSimul) =  p(ii, iSimul) - p(ii-1, iSimul);
        end
    elseif strcmpi(policyTargetRule, 'inflation')
        pii(ii, iSimul) =  sum(nrGap(ii:ii+nFuture));
        if ii == 1
            p(ii, iSimul) =  pii(ii, iSimul) - 0;
        else
            p(ii, iSimul) =  p(ii-1, iSimul) + pii(ii, iSimul);
        end
    end
    
end

%% Plot graphs

nPlotPeriods = 20;

f = figure;

previousInterpreter = pub_GraphSetInterpreter('Latex');

if strcmpi(policyTargetRule, 'priceLevel')
    yData = p;
    ruleLabel = '$p_t$';
    simGraphName = ['graph_WicksellianTheoryRiskyNatural_', policyTargetRule, '_phi', replace(num2str(phi), '.', ''), '.png'];
elseif strcmpi(policyTargetRule, 'inflation')
    yData = pii;
    ruleLabel = '$\pi_t$';
    simGraphName = ['graph_WicksellianTheoryRiskyNatural_', policyTargetRule, '_phi', replace(num2str(phiPii), '.', ''), '.png'];
end

tl = tiledlayout(3,2, 'tilespacing', 'compact', 'padding', 'compact');

iPlot = 1;
nexttile();
graph = plot(100 * yData(:,1), 'LineWidth', 3);
graph(1).LineStyle = '-';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("Steady state", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-30,30]);
yticks([-30:15:30]);
pub_GraphDrawZeroAxis(graph);
grid 'on';

iPlot = 2;
nexttile();
graph = plot(100*[yData(:,2), yData(:,3)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$\mathcal{D}^r$ goes up from 10\% to 15\% from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('');
xlim([0,nPlotPeriods]);
ylim([-30,30]);
yticks([-30:15:30]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated shocks', 'MIT shocks', 'Interpreter', 'latex');
grid 'on';

iPlot = 3;
nexttile();
graph = plot(100*[yData(:,4), yData(:,5)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^{Good}$ goes up by 50\% from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-30,30]);
yticks([-30:15:30]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated shocks', 'MIT shocks', 'Interpreter', 'latex');
grid 'on';

iPlot = 4;
nexttile();
graph = plot([100*yData(:,6), yData(:,7)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^{Bad}$ goes up from -60\% to -50\% from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-30,30]);
yticks([-30:15:30]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated shocks', 'MIT shocks', 'Interpreter', 'latex');
grid 'on';

iPlot = 5;
nexttile();
graph = plot(100*[yData(:,8), yData(:,9)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$\bar{\iota}$ goes up by 25\% from t = 1 to 5", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-30,30]);
yticks([-30:15:30]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated shocks', 'MIT shocks', 'Interpreter', 'latex');
grid 'on';

iPlot = 6;
nexttile();
graph = plot(100*[yData(:,10), yData(:,11)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^n$ is a stochastic process", 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel(ruleLabel, 'Interpreter', 'latex');
xlabel('', 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-30,30]);
yticks([-30:15:30]);
pub_GraphDrawZeroAxis(graph);
legend('fixed intercept (ignoring risk)', '$r^n$ in the intercept', 'Interpreter', 'latex');
grid 'on';

%edu_Suptitle(['Neo-Wicksellian theory of the price level with risky natural ($\phi$ = ', num2str(phi), ')'], ...
%    'FontSize', 20, 'Interpreter', 'latex');

pub_GraphSetInterpreter(previousInterpreter);

set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);

