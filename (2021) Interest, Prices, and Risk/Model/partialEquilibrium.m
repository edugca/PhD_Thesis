%% Wicksellian model

isFinal = false;

if isFinal
    imagesFolder = '/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Images/Partial Equilibrium/';
else
    imagesFolder = 'Model/Images/';
end

%%% To replicate the Paper's figures, change phi (or phiPii) and set the respective policyTarget rule
phi         = 0.2; % price-level targeting: monetary policy parameter > 0
phiPii      = 1.2; % inflation targeting:monetary policy parameter > 1

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

fontSize        = 14;
titleFontSize   = 16;

%% Partial equilibrium without policy default probability

%%% Steady state
iSimul       = 1;
rN(:)        = 0;
iota(:)      = 0;
defPol(:)    = 0;

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

%%% Jump in rN
iSimul        = iSimul + 1;
rN(:)         = 0;
iota(:)       = 0;
defPol(:)    = 0;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            rN(j)        = log(1.2);
        end
        
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

%%% Jump in rN (MIT shocks)
iSimul        = iSimul + 1;
rN(:)         = 0;
iota(:)       = 0;
defPol(:)    = 0;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            rN(ii)        = log(1.2);
        end
        
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
rN(:)        = 0;
iota(:)      = 0;
defPol(:)    = 0;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if j >= 1 && j <= 5
            iota(j)        = log(1.2);
        end
        
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
rN(:)        = 0;
iota(:)      = 0;
defPol(:)    = 0;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        % Shock rule
        if ii >= 1 && ii <= 5 && j == ii
            iota(ii)        = log(1.2);
        end
        
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
rN(:)        = real(log(1 + 0.02*randn(length(rN),1)));
iota(:)      = 0;
defPol(:)    = 0;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(ii) - iota(ii));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(ii) - iota(ii));
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
defPol(:)    = 0;

nrGap(:) = 0;
for ii=1:nPeriods
    for j=ii:(ii+nFuture)
        if strcmpi(policyTargetRule, 'priceLevel')
            nrGap(j)   = alpha^(j-ii+1) * (rN(ii) - iota(ii));
        elseif strcmpi(policyTargetRule, 'inflation')
            nrGap(j)   = alphaPii^(j-ii+1) * (rN(ii) - iota(ii));
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
    simGraphName = ['graph_WicksellianTheory_', policyTargetRule, '_phi', replace(num2str(phi), '.', ''), '.png'];
elseif strcmpi(policyTargetRule, 'inflation')
    yData = pii;
    ruleLabel = '$\pi_t$';
    simGraphName = ['graph_WicksellianTheory_', policyTargetRule, '_phi', replace(num2str(phiPii), '.', ''), '.png'];
end

tl = tiledlayout(2,2, 'tilespacing', 'compact', 'padding', 'compact');

iPlot = 1;
nexttile();
graph = plot(100 * yData(:,1), 'LineWidth', 3);
graph(1).LineStyle = '-';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("Steady state", 'FontSize', titleFontSize, 'Interpreter', 'latex')
ylabel(ruleLabel, 'FontSize', fontSize, 'Interpreter', 'latex');
xlabel('Time', 'FontSize', fontSize, 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-100,100]);
pub_GraphDrawZeroAxis(graph);
grid 'on';

iPlot = 2;
nexttile();
graph = plot(100*[yData(:,2), yData(:,3)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^n$ goes up by 20\% from t = 1 to 5", 'FontSize', titleFontSize, 'Interpreter', 'latex');
ylabel(ruleLabel, 'FontSize', fontSize, 'Interpreter', 'latex');
xlabel('Time', 'FontSize', fontSize, 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-100,100]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated shocks', 'MIT shocks', 'FontSize', fontSize, 'Interpreter', 'latex');
grid 'on';

iPlot = 3;
nexttile();
graph = plot(100*[yData(:,4), yData(:,5)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title('$\bar{\iota}$ goes up by 20\% from t = 1 to 5', 'FontSize', titleFontSize, 'Interpreter', 'latex');
ylabel(ruleLabel, 'FontSize', fontSize, 'Interpreter', 'latex');
xlabel('Time', 'FontSize', fontSize, 'Interpreter', 'latex');
xlim([0,nPlotPeriods]);
ylim([-100,100]);
pub_GraphDrawZeroAxis(graph);
legend('anticipated shocks', 'MIT shocks', 'FontSize', fontSize, 'Interpreter', 'latex');
grid 'on';

iPlot = 4;
nexttile();
graph = plot(100*[yData(:,6), yData(:,7)], 'LineWidth', 3);
graph(1).LineStyle = '-';
graph(2).LineStyle = ':';
set(gca,'FontSize',fontSize); % Scale fontsize of axes
title("$r^n$ is a stochastic process", 'FontSize', titleFontSize, 'Interpreter', 'latex');
ylabel(ruleLabel, 'FontSize', fontSize, 'Interpreter', 'latex');
xlabel('Time', 'FontSize', fontSize);
xlim([0,nPlotPeriods]);
ylim([-100,100]);
pub_GraphDrawZeroAxis(graph);
legend('fixed intercept', '$r^n$ in the intercept', 'FontSize', fontSize, 'Interpreter', 'latex');
grid 'on';

%edu_Suptitle(['Neo-Wicksellian theory of the price level ($\phi$ = ', num2str(phi), ')'], ...
%                'FontSize', 20, 'Interpreter', 'latex');

pub_GraphSetInterpreter(previousInterpreter);

set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
saveas(f,[imagesFolder, simGraphName]);