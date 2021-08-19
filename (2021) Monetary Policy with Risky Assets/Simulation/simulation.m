%%% Paper: Monetary Policy with Risky Assets
%%%
%%% Author: Eduardo G. C. Amaral
%%% Creation date: November 9, 2020
%%% Last modification: November 9, 2020
%%%
%%% Wicksellian model: simulation of correlation
%%%
%%% In case you find any mistake in this code, do not hesitate in sending
%%% me the heads up :)
%%%

%% Housekeeping

clear all;
clc;

isFinal = false;

if isFinal
    imagesFolder = '/Users/Eduardo/Dropbox/PUC-Rio/Tese/Natural Interest Rates/Paper 1/Images/Partial Equilibrium/';
else
    imagesFolder = '/Users/Eduardo/OneDrive/MATLAB/My Thesis/Paper 1/Wicksellian Model/Images/';
end

%% Set up

%%%% Change the policyTargetRule to replicate the Paper's figures

policyTargetRule = 'inflation'; % priceLevel, inflation

%%%%%%%%%%%%%%%%%%%%%%%

nPeriods        = 1;
nFuture         = 100;
nSimulations    = 10000;
p               = nan(nPeriods, nSimulations);
pii             = nan(nPeriods, nSimulations);
rN              = nan(nPeriods + nFuture, nSimulations);
iota            = nan(nPeriods + nFuture, nSimulations);
nrGap           = nan(nPeriods + nFuture, nSimulations);
defPol          = nan(nPeriods + nFuture, nSimulations);
nrBad           = nan(nPeriods + nFuture, nSimulations);

%% Run simulations

simData = struct();

% Initiate the random seeder
rng(1900);

% List of set-ups
if strcmp(policyTargetRule, 'priceLevel')
    phiList     = {0.2, 0.5, 1.0};
elseif strcmp(policyTargetRule, 'inflation')
    phiList     = {1.2, 1.5, 2.0}; 
end
correlList  = {-1, -0.5, 0, 0.5 ,1};

for iCorrel = 1:length(correlList)
    
    p_Vec       = containers.Map();
    pii_Vec     = containers.Map();
    rN_Vec      = containers.Map();
    defPol_Vec  = containers.Map();
    
    for iPhi = 1:length(phiList)
    
        iExp = 0;
        expNames = {};
        correl_rN_defPol = correlList{iCorrel};
        
        if strcmp(policyTargetRule, 'priceLevel')
            phi     = phiList{iPhi};    % price-level targeting: monetary policy parameter > 0
            alpha       = 1/(1+phi);
        elseif strcmp(policyTargetRule, 'inflation')
            phiPii  = phiList{iPhi};    % inflation targeting:monetary policy parameter > 1
            alphaPii    = 1/phiPii;
        end
        
        if strcmp(policyTargetRule, 'priceLevel')
            thisPhi = phi;
        elseif strcmp(policyTargetRule, 'inflation')
            thisPhi = phiPii;
        end
        
        p_Vec(num2str(thisPhi))        = [];
        pii_Vec(num2str(thisPhi))      = [];
        rN_Vec(num2str(thisPhi))       = [];
        defPol_Vec(num2str(thisPhi))   = [];
        
        simulation_experiences;
    end
    
    simData.(['Correl_' num2str(iCorrel)]).p_Vec = p_Vec;
    simData.(['Correl_' num2str(iCorrel)]).pii_Vec = pii_Vec;
    simData.(['Correl_' num2str(iCorrel)]).rN_Vec = rN_Vec;
    simData.(['Correl_' num2str(iCorrel)]).defPol_Vec = defPol_Vec;
    
end

%% Plot graphs

% Phis to plot
phisToPlot = cellfun(@(x) num2str(x), phiList, 'UniformOutput', false);

combList        = {    'p_Vec',     'defPol_Vec';
                       'pii_Vec',   'defPol_Vec'
                  };
        
combNamesList   = {    '$\hat{p}_t$',     '$E_t \mathcal{D}^{Policy}_{t+1}$';
                       '$\pi_t$',   '$E_t \mathcal{D}^{Policy}_{t+1}$'
                  };
graphNamesList  = {    'P_defPol';
                       'pii_defPol'
                  };
              
fontSize = 12;
              
% Iterate variable combinations
for iComb = 1:size(combList, 1)
    % Iterate correlations
    for iCorrel = 1:length(correlList)
        f = figure;
        previousInterpreter = pub_GraphSetInterpreter('Latex');
        nLins = length(phisToPlot);
        nCols = iExp + 1;

        % [left right top bottom]  1 = 100%
        offsetSpecs = [0.1, 0.05 , 0.05 , 0.1];
        % [width height] 1 = 100%
        scaleSpecs = [0.9, 0.9];
        pos = iosr.figures.subfigrid(nLins, nCols, offsetSpecs, scaleSpecs);

        markerSize = 5;
        for iPhi = 1:length(phisToPlot)

            plotX = simData.(['Correl_' num2str(iCorrel)]).rN_Vec(phisToPlot{iPhi});
            plotY = simData.(['Correl_' num2str(iCorrel)]).defPol_Vec(phisToPlot{iPhi});

            dataX = plotX(:, 1);
            dataY = plotY(:, 1);

            % Plot the correlation between shocks
            sp = subplot('position',pos(iPhi,:,1));
            %subplot(nLins, nCols, 1 + (iPhi-1)*nCols);
            graph = scatter(dataX, dataY, markerSize, 'filled');
            set(sp,'FontSize', fontSize); % Scale fontsize of axes
            
            pCoeff  = polyfit(dataX, dataY,1); % This will generate the coefficients of polynomial of degree(1)
            yFitted = polyval(pCoeff, dataX); % This will give the fitted  values for the desired values.
            hold on;
            plot(dataX, yFitted,'k-');
            legend('boxoff');
            
            % Calculate correlations (H0 = There is no correlation; H1 = correlation is different from 0)
            [corrPearson, pPearson]   =  corr(dataX, dataY, 'Rows', 'complete', 'Type', 'Pearson', 'Tail', 'both');
            [corrSpearman, pSpearman] =  corr(dataX, dataY, 'Rows', 'complete', 'Type', 'Spearman', 'Tail', 'both');
            leg = legend({'', [num2str(corrPearson,'%.2f'), ' (', char(num2str(pPearson, '%.2f')), ')']}, 'location', 'best');

            if iPhi == 1
                title('$Corr \left( r^n_t, E_t\mathcal{D}^{Policy}_{t+1} \right)$', 'Interpreter', 'latex');
            end
            if strcmp(policyTargetRule, 'priceLevel')
                ylabel({['$\phi = ' , phisToPlot{iPhi}, '$']; combNamesList{iComb, 2}});
            elseif strcmp(policyTargetRule, 'inflation')
                ylabel({['$\phi^{\pi} = ' , phisToPlot{iPhi}, '$']; combNamesList{iComb, 2}});
            end
            if iPhi == length(phisToPlot)
                xlabel('$r^n_t$');
            end
            %xlim([0,nPlotPeriods]);
            %ylim([-25,25]);

            for iPlot = 1:iExp

                plotX = simData.(['Correl_' num2str(iCorrel)]).( combList{iComb,1} )( phisToPlot{iPhi} );
                plotY = simData.(['Correl_' num2str(iCorrel)]).( combList{iComb,2} )( phisToPlot{iPhi} );

                dataX = plotX(:, iPlot);
                dataY = plotY(:, iPlot);

                sp = subplot('position',pos(iPhi,:,1 + iPlot));
                %subplot(nLins, nCols, 1 + iPlot + (iPhi-1)*nCols);
                set(sp,'FontSize', fontSize); % Scale fontsize of axes
                graph = scatter(dataX, dataY, markerSize, 'filled');
                set(sp,'FontSize', fontSize); % Scale fontsize of axes
                
                pCoeff  = polyfit(dataX, dataY,1); % This will generate the coefficients of polynomial of degree(1)
                yFitted = polyval(pCoeff, dataX); % This will give the fitted  values for the desired values.
                hold on;
                plot(dataX, yFitted,'k-');
                legend('boxoff');
                
                % Calculate correlations (H0 = There is no correlation; H1 = correlation is greater than 0)
                [corrPearson, pPearson]   =  corr(dataX, dataY, 'Rows', 'complete', 'Type', 'Pearson', 'Tail', 'right');
                [corrSpearman, pSpearman] =  corr(dataX, dataY, 'Rows', 'complete', 'Type', 'Spearman', 'Tail', 'right');
                leg = legend({'', [num2str(corrPearson,'%.2f'), ' (', char(num2str(pPearson, '%.2f')), ')']}, 'location', 'best');

                % Remove tick labels
                yticklabels({});

                if iPhi == 1
                    title(expNames{iPlot}, 'Interpreter', 'latex');
                end
                %ylabel('$E_t \mathcal{D}_{t+1}$');
                if iPhi == length(phisToPlot)
                    xlabel(combNamesList{iComb, 1});
                end
                %xlim([0,nPlotPeriods]);
                %ylim([-25,25]);
                %set(gca,'FontSize',16); % Scale fontsize of axes
                
            end
            
            %set(gca,'FontSize',16); % Scale fontsize of axes
        end

        % Link x axes of each column
        for iCol = 1:nCols
            listAxes = [];
            for iLin = 1:nLins
                %listAxes = [listAxes, subplot(nLins, nCols, iAx)];
                listAxes = [listAxes, subplot('position',pos(iLin,:,iCol))];
            end
            linkaxes(listAxes, 'x');
        end

        % Link y axes of each line
        for iLin = 1:nLins
            listAxes = [];
            for iCol = 1:nCols
                %listAxes = [listAxes, subplot(nLins, nCols, iAx)];
                listAxes = [listAxes, subplot('position',pos(iLin,:,iCol))];
            end
            linkaxes(listAxes, 'y');
        end

        %edu_Suptitle(['Neo-Wicksellian theory of the price level with policy asset default ($\phi$ = ', num2str(phi), ')'], ...
        %    'FontSize', 20, 'Interpreter', 'latex');

        simGraphName = ['graph_WicksellianTheoryDefault_rule' policyTargetRule '_' graphNamesList{iComb} '_Correlation_byPhi_correl', replace(num2str(correlList{iCorrel}), '.', ''), '.png'];
        set(gcf, 'Position',  [100, 100, 1000, 600]); % resize figure
        saveas(f,[imagesFolder, simGraphName]);

        pub_GraphSetInterpreter(previousInterpreter);
    end
end