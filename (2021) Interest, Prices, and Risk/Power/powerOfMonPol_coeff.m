%%% Paper: Monetary Policy with Risky Assets
%%%
%%% Author: Eduardo G. C. Amaral
%%% Creation date: November 9, 2020
%%% Last modification: November 9, 2020
%%%
%%% Wicksellian model: monetary policy power
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
    imagesFolder = 'Power/Images/';
end

%% Parameterization

% Pick the regime
regime = 'Inflation'; % 'PriceLevel' or 'Inflation'

if strcmp(regime, 'PriceLevel')
    phiRegime = [0.1, 0.5, 1.0];
elseif strcmp(regime, 'Inflation')
    phiRegime = [1.1, 1.5, 2.0];
end

%% Function

nSims   = 100;
nVec    = 40;
nPlot   = 40;

maxDef  = 0.20;
maxDef  = maxDef * 1000000;

for corrCoeff = [-1, -0.5, 0, 0.5, 1]
    
    f = figure;
    previousInterpreter = pub_GraphSetInterpreter('Latex');
    tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
    
    mu      = [0.10 0.04]; % Mean of epected default probability and rN
    Sigma   =   [
                    0.01.^2                             corrCoeff*(0.01)*(0.005); 
                    corrCoeff*(0.01)*(0.005)            0.005.^2
                ];
    
    for iSim = 1:nSims
       R(:,:,iSim) = abs(mvnrnd(mu, Sigma, nVec)); 
    end
    
    iPhi = 0;
    for phi = phiRegime
        
        iPhi = iPhi + 1;

        for iSim = 1:nSims

            v_def   = R(:,1,iSim) ;
            rn      = R(:,2,iSim);
            iota    = 0.04 .* ones(nVec,1);
            
            if strcmp(regime, 'PriceLevel')
                [directRiskFree, gammaRiskFree, directRisky, gammaRisky] = fPowerMonPol_PLTarg(phi, v_def);
            elseif strcmp(regime, 'Inflation')
                [directRiskFree, gammaRiskFree, directRisky, gammaRisky] = fPowerMonPol_InfTarg(phi, v_def);
            end
            v_riskFree  = directRiskFree .* gammaRiskFree;
            v_Risky     = directRisky .* gammaRisky;

            nexttile((iPhi-1)*4 + 1);
            hold on;
            p = scatter(rn(1:nPlot), v_def(1:nPlot), 0.05);
            p.LineWidth = 0.01;
            p.MarkerEdgeColor = 'b';
            p.MarkerFaceColor = [0 0.5 0.5];
    %         p(1).LineWidth = 0.01;
    %         p(1).Color     = '#d8dcd6';
            if iSim == nSims
                pub_GraphDrawZeroAxis(p);
                set(gca,'FontSize',12); % Scale fontsize of axes
                if strcmp(regime, 'PriceLevel')
                    ylabel({['$\phi = ' , num2str(phi), '$']; '$E_t \mathcal{D}_{t+1}\delta_{t+1}$'}, 'Interpreter', 'latex');
                elseif strcmp(regime, 'Inflation')
                    ylabel({['$\phi^{\pi} = ' , num2str(phi), '$']; '$E_t \mathcal{D}_{t+1}\delta_{t+1}$'}, 'Interpreter', 'latex');
                end
                
                if iPhi == 1
                    title("Scatter plot", 'Interpreter', 'latex')
                elseif iPhi == 3
                   xlabel('$r^n_t$', 'Interpreter', 'latex'); 
                end
            end

            nexttile((iPhi-1)*4 + 2);
            hold on;
            p = plot(1:nPlot, gammaRisky(1:nPlot) - gammaRiskFree(1:nPlot));
            p(1).LineWidth = 0.01;
            p(1).Color     = '#d8dcd6';
            if iSim == nSims
                pub_GraphDrawZeroAxis(p);
                %ylim([-0.5 1]);
                set(gca,'FontSize',12); % Scale fontsize of axes
                ylabel("Difference", 'Interpreter', 'latex');
                
                if iPhi == 1
                    title("$\Upsilon^{Risky}_{t,j+1} - \Upsilon^{RF}_{t,j+1}$", 'Interpreter', 'latex')
                elseif iPhi == 3
                    xlabel('term period', 'Interpreter', 'latex');
                end
            end

            nexttile((iPhi-1)*4 + 3);
            hold on;
            p = plot(1:nPlot, v_Risky(1:nPlot) - v_riskFree(1:nPlot));
            p(1).LineWidth = 0.01;
            p(1).Color     = '#d8dcd6';
            if iSim == nSims
                pub_GraphDrawZeroAxis(p);
                %ylim([-0.5 1]);
                set(gca,'FontSize',12); % Scale fontsize of axes
                ylabel("Difference", 'Interpreter', 'latex');
                
                if iPhi == 1
                    title("$\Upsilon^{Risky}_{t,j+1}\left(1 - E_t \mathcal{D}_{t+j+1}\delta_{t+j+1}\right) - \Upsilon^{RF}_{t,j+1}$", 'Interpreter', 'latex')
                elseif iPhi == 3
                    xlabel('term period', 'Interpreter', 'latex');
                end
            end

            nexttile((iPhi-1)*4 + 4);
            hold on;
            p = plot(1:nPlot, 100 .* ( gammaRisky(1:nPlot).*(rn(1:nPlot) - directRisky(1:nPlot).*iota(1:nPlot)) - gammaRiskFree(1:nPlot).*(rn(1:nPlot) - directRiskFree(1:nPlot).*iota(1:nPlot)) ));
            p(1).LineWidth = 0.01;
            p(1).Color     = '#d8dcd6';
            if iSim == nSims
                pub_GraphDrawZeroAxis(p);
                set(gca,'FontSize',12); % Scale fontsize of axes
                ylabel("%", 'Interpreter', 'latex');
               
                if iPhi == 1
                    if strcmp(regime, 'PriceLevel')
                        title("$\Delta p_t$ per term", 'Interpreter', 'latex')
                    elseif strcmp(regime, 'Inflation')
                        title("$\Delta \pi_t$ per term", 'Interpreter', 'latex')
                    end
                elseif iPhi == 3
                    xlabel('term period', 'Interpreter', 'latex');
                end
            end
        end
    end
    
    linkaxes([nexttile(1) nexttile(5) nexttile(9)]);
    linkaxes([nexttile(2) nexttile(6) nexttile(10)]);
    linkaxes([nexttile(3) nexttile(7) nexttile(11)]);
    linkaxes([nexttile(4) nexttile(8) nexttile(12)]);
    
    pub_GraphSetInterpreter(previousInterpreter);

    set(gcf, 'Position',  [100, 100, 1000, 800]); % resize figure
    simGraphName = ['monPolPower_simulation_' regime '_corr_' num2str(corrCoeff) '.png'];
    saveas(f,[imagesFolder, simGraphName]);

end

%% Function

% Price-level targeting
function [directRiskFree, gammaRiskFree, directRisky, gammaRisky] = fPowerMonPol_PLTarg(phi, def)

gammaRiskFree   = ones(length(def), 1);
directRiskFree  = ones(length(def), 1);
gammaRisky      = ones(length(def), 1);
directRisky     = ones(length(def), 1);

for t = 1:length(def)
 
    if t ~= 1
        gammaRiskFree(t)    = gammaRiskFree(t-1)* ((1 + phi)^(-1)) ;
        gammaRisky(t)       = (1-def(t)) * gammaRisky(t-1)* ((1 + (1 - def(t))*phi)^(-1)) ;
    else
        gammaRiskFree(t)    = ((1 + phi)^(-1)) ;
        gammaRisky(t)       = ((1 + (1 - def(t))*phi)^(-1)) ;
    end
    
    directRiskFree(t)    = 1 ;
    directRisky(t)       = (1-def(t)) ;
        
end

riskFree    = directRiskFree .* gammaRiskFree;
risky       = directRisky .* gammaRisky;

end

% Inflation targeting
function [directRiskFree, gammaRiskFree, directRisky, gammaRisky] = fPowerMonPol_InfTarg(phi, def)

gammaRiskFree   = ones(length(def), 1);
directRiskFree  = ones(length(def), 1);
gammaRisky      = ones(length(def), 1);
directRisky     = ones(length(def), 1);

for t = 1:length(def)
 
    if t ~= 1
        gammaRiskFree(t)    = gammaRiskFree(t-1)* ((phi)^(-1)) ;
        gammaRisky(t)       = (1-def(t)) * gammaRisky(t-1)* (((1 - def(t))*phi)^(-1)) ;
    else
        gammaRiskFree(t)    = ((phi)^(-1)) ;
        gammaRisky(t)       = (((1 - def(t))*phi)^(-1)) ;
    end
    
    directRiskFree(t)    = 1 ;
    directRisky(t)       = (1-def(t)) ;
        
end

riskFree    = directRiskFree .* gammaRiskFree;
risky       = directRisky .* gammaRisky;

end