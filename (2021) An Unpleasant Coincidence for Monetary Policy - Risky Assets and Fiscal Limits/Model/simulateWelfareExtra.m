%% Plot welfare results

close all;

simAllMdls_Plot = cellfun(@(x) x(iParams,jParams), simAllMdls, 'UniformOutput', false);
simAllMdls_R_Plot = cellfun(@(x) x(iParams,jParams), simAllMdls_R, 'UniformOutput', false);
simAllMdls_NR_Plot = cellfun(@(x) x(iParams,jParams), simAllMdls_NR, 'UniformOutput', false);

% Create color map
simNormalized =  normalize(cell2mat([simAllMdls_Plot; simAllMdls_R_Plot; simAllMdls_NR_Plot]), 'range', [0 1]);
colorMapLims = [floor(min(min(simNormalized))) ceil(max(max(simNormalized)))];

% create a default color map ranging from red to light pink
colorLength = 100;
color1    = [1, 0, 0];
color2    = [1, 1, 1];
myColorMap = [linspace(color1(1),color2(1),colorLength)', linspace(color1(2),color2(2),colorLength)', linspace(color1(3),color2(3),colorLength)'];
color1    = [1, 1, 1];
color2    = [0, 0, 1];
myColorMap = [myColorMap; [linspace(color1(1),color2(1),colorLength)', linspace(color1(2),color2(2),colorLength)', linspace(color1(3),color2(3),colorLength)']];


% cLim = [min(-vec([simAllMdls{:}])), max(-vec([simAllMdls{:}]))];
% centerPoint = prctile(-vec([simAllMdls{:}]), 0.5);
% scalingIntensity = 5;
% inc = 100;
% myColors = [1, 0, 0; 1, 1, 1; 0, 0, 1];
% [myColorMap, ticksMap, tickLabelsMap] = eduColorMapNonlinear(myColors, centerPoint, cLim, scalingIntensity, inc);

if iscell(paramNames{1})
    str_paramNames_1 = paramNames{1}{1};
else
    str_paramNames_1 = paramNames{1};
end

if iscell(paramNames{2})
    str_paramNames_2 = paramNames{2}{1};
else
    str_paramNames_2 = paramNames{2};
end


for iMdl = 1:length(mdlSimVector)
    config  = riseConfigVector(iMdl);
    modelName   = riseModelNames{iMdl};
    defProcess  = riseConfigVector(iMdl).conf_defRProcessType;
    
    %simResults   = (-simAllMdls{iMdl}) - (-simAllMdls{refMdl});
    %simResultsR  = (-simAllMdls_R{iMdl}) - (-simAllMdls{refMdl});
    %simResultsNR = (-simAllMdls_NR{iMdl}) - (-simAllMdls{refMdl});
    
    simResults   = (-simAllMdls{iMdl});
    simResultsR  = (-simAllMdls_R{iMdl});
    simResultsNR = (-simAllMdls_NR{iMdl});
    
    %nLins = size(simAllMdls_Plot{iMdl},1);
    %nCols = size(simAllMdls_Plot{iMdl},2);
    %simResults   = simNormalized(1:nLins, (nCols*(iMdl-1) + 1):(nCols*iMdl));
    %simResultsR  = simNormalized((nLins+1):(2*nLins),  (nCols*(iMdl-1) + 1):(nCols*iMdl));
    %simResultsNR = simNormalized((2*nLins+1):(3*nLins), (nCols*(iMdl-1) + 1):(nCols*iMdl));

    if length(size(simResults)) == 2
        %%%% Aggregate welfare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        previousInterpreter = pub_GraphSetInterpreter('latex');
        h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResults));
        colormap(myColorMap);
        caxis(colorMapLims);
        xlabel(paramPlotNames{2});
        ylabel(paramPlotNames{1});
        set(gca, 'FontSize', 12);
        
        yLims = ylim;
        ylim([yLims(1) 1]);
        
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['Aggregate_withGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        else
                simGraphName = ['Aggregate_withoutGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        end
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%% Ricardian welfare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        previousInterpreter = pub_GraphSetInterpreter('latex');
        h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResultsR));
        colormap(myColorMap);
        caxis(colorMapLims);
        xlabel(paramPlotNames{2});
        ylabel(paramPlotNames{1});
        set(gca, 'FontSize', 12);
        
        yLims = ylim;
        ylim([yLims(1) 1]);
        
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['Ricardian_withGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        else
                simGraphName = ['Ricardian_withoutGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        end
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%% Non-Ricardian welfare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        previousInterpreter = pub_GraphSetInterpreter('latex');
        h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResultsNR));
        colormap(myColorMap);
        caxis(colorMapLims);
        xlabel(paramPlotNames{2});
        ylabel(paramPlotNames{1});
        set(gca, 'FontSize', 12);
        
        yLims = ylim;
        ylim([yLims(1) 1]);
        
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['NonRicardian_withGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        else
                simGraphName = ['NonRicardian_withoutGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        end
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
        %%%% Line Plot welfare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        nCol = 3;
        nLin = ceil(length(paramValues{2}) / nCol);
        tl = tiledlayout(nLin,nCol);
        tl.Padding = 'compact';
        tl.TileSpacing = 'compact';
        previousInterpreter = pub_GraphSetInterpreter('latex');
        
        lab = paramValues{1};
        jParams = 1:length(paramValues{2});
        
        for iParam = 1:length(paramValues{2})
            h1 = nexttile();
            h = plot(paramValues{2}(jParams), simResults(iParam,jParams), 'DisplayName', 'Aggregate');
            h.LineStyle = '-';
            h.LineWidth = 2;
            hold on;
            h = plot(paramValues{2}(jParams), simResultsR(iParam,jParams), 'DisplayName', 'Ricardian');
            h.LineStyle = '--';
            h.LineWidth = 2;
            hold on;
            h = plot(paramValues{2}(jParams), simResultsNR(iParam,jParams), 'DisplayName', 'Non-Ricardian');
            h.LineStyle = '-.';
            h.LineWidth = 2;
            grid('on');
            xlabel(paramTexPlotNames{2});
            ylabel('Welfare');
            title(['$\phi = $' num2str(lab(iParam))]);
            set(gca, 'FontSize', 12);
            xlim([paramValues{2}(jParams(1)) paramValues{2}(jParams(end))]);
            
            xl = xline(simStruct.params.(paramNames{2}), 'DisplayName', 'Calibration');
            xl.LineStyle = ':';
            xl.LineWidth = 2;
            xl.Color = 'black';
            
            [maxVal,idx] = max(simResults(iParam,jParams));
            hold on
            yMax = plot(paramValues{2}(jParams(idx)), maxVal,'r*', 'DisplayName', 'Maximum');
            yMax.LineWidth  = 2;
            yMax.Color      = 'red';
        end
        lg = legend('Location', 'best', 'Orientation', 'horizontal', 'FontSize', 12);
        lg.Layout.Tile = 'North';

        linkaxes(tl.Children, 'y');
        
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['Welfare_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '_Selected' '.png'];
        else
                simGraphName = ['Welfare_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '_Selected' '.png'];
        end
        set(gcf, 'Position',  [100, 100, 1000, 700]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%% TABLE
        
%         t = table();
%         t.Aggregate     = normalize(simResults, 'scale');
%         t.Ricardian     = normalize(simResultsR, 'scale');
%         t.NonRicardian  = normalize(simResultsNR, 'scale');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif length(size(simResults)) == 3
        
        %%%% Difference in welfare between NR and R agents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        iParam = 1;
        previousInterpreter = pub_GraphSetInterpreter('latex');
        h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResultsNR(:,:,iParam) - simResultsR(:,:,iParam)));
        colormap(myColorMap);
        caxis(colorMapLims);
        xlabel(paramPlotNames{2});
        ylabel(paramPlotNames{1});
        set(gca, 'FontSize', 12);
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['Aggregate_withGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        else
                simGraphName = ['Aggregate_withoutGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        end
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%% Difference in welfare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        iParam = 4;
        previousInterpreter = pub_GraphSetInterpreter('latex');
        h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResultsNR(:,:,iParam) - simResultsR(:,:,iParam)));
        colormap(myColorMap);
        caxis(colorMapLims);
        xlabel(paramPlotNames{2});
        ylabel(paramPlotNames{1});
        set(gca, 'FontSize', 12);
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['Aggregate_withGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        else
                simGraphName = ['Aggregate_withoutGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        end
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%% Difference in welfare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = figure;
        iParam = 7;
        previousInterpreter = pub_GraphSetInterpreter('latex');
        h = heatmap(paramValues{2}, fliplr(paramValues{1}), flipud(simResultsNR(:,:,iParam) - simResultsR(:,:,iParam)));
        colormap(myColorMap);
        caxis(colorMapLims);
        xlabel(paramPlotNames{2});
        ylabel(paramPlotNames{1});
        set(gca, 'FontSize', 12);
        pub_GraphSetInterpreter(previousInterpreter);
        if config.conf_hasGovernment
                simGraphName = ['Aggregate_withGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        else
                simGraphName = ['Aggregate_withoutGov_paramStability_' str_paramNames_1 '_and_' str_paramNames_2 '_' config.conf_policyRule '.png'];
        end
        set(gcf, 'Position',  [100, 100, 800, 800]); % resize figure
        exportgraphics(f, ...
        strjoin({pathImages 'Welfare' title_stickyPrices simGraphName}, filesep));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
%     % Add benchmark parameterization
%     ax = axes;
%     scatter(ax, paramRef{2}, paramRef{1},'LineWidth',2)
%     ax.Color = 'none'; % make invisible the background of the scatter plot

end


%% Save results

simGraphName = ['Welfare ' char(replace( string(datetime), ':', '-')) '.mat'];
filePath = strjoin({pathSaved 'Welfare' title_stickyPrices simGraphName}, filesep);

save(filePath);