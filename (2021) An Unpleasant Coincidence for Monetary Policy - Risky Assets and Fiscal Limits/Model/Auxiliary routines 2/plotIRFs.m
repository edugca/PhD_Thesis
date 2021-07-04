%% Plot IRFs
function graph = plotIRFs(figureRef, plotVars, periods, suptitleStr, ...
    legendStr, irfLevel, irfLevelChg, irfPercChg, nPlotCols, ...
    plotVarNames, legendInSubPlot, setFirstPeriodToNaN, customizeLines)
    
    %%% PARAMETERS
    %
    % customizeLines = struct() => user-defined linestyle; [] this function linestyle; false => Matlab linestyle
    %
    %%%%%%%%%%%%%%%%%%
    
    
    % minimum % change to not be treated as ZERO
    tol = 10e-15;
    
    nPlotVars = length(plotVars);
    auxDesloc = legendInSubPlot;
    if isempty(nPlotCols)
        nPlotLines = 2;
        nPlotCols = ceil((nPlotVars + auxDesloc) / nPlotLines);
    else
        nPlotLines = ceil((nPlotVars + auxDesloc) / nPlotCols);
    end
    
    % Define figure object
    if isempty(figureRef)
        graph = figure;
        
        tl = tiledlayout(nPlotLines, nPlotCols);
        tl.TileSpacing  = 'compact';
        tl.Padding      = 'compact';
        
        % Check whether legend in first subplot
        desloc = 0;
        if legendInSubPlot > 0
            desloc = legendInSubPlot;
            subGraph_legend = nexttile(1);
        end
        
        % Create subGraphs and attribute a tag to each
        for i=1:nPlotVars
            subGraph        = nexttile(i + desloc);
            subGraph.Tag    = plotVars{i};
        end
    elseif isa(figureRef, 'matlab.graphics.axis.Axes') && length(figureRef) > 1
        
        graph = figureRef;
        
        % Check whether legend in first subplot
        desloc = 0;
        if legendInSubPlot > 0
            desloc = legendInSubPlot;
            subGraph_legend = graph(1);
        end
        
        % Create subGraphs and attribute a tag to each
        for i=1:nPlotVars
            subGraph = graph(i + desloc);
            subGraph.Tag = plotVars{i};
        end
    else
        graph = figureRef;
    end
    
    for i=1:nPlotVars
        
        % Check whether graph is axes or figure
        if isa(graph,'matlab.graphics.axis.Axes') && length(figureRef) == 1
            if isa(graph.Parent, 'matlab.graphics.layout.TiledChartLayout')
                set(graph.Parent.Parent,'CurrentAxes',figureRef);
            else
                set(graph.Parent,'CurrentAxes',figureRef);
            end
        elseif isa(graph,'matlab.graphics.axis.Axes') && length(figureRef) > 1
            subGraph = graph.findobj('Tag', plotVars{i});
            set(subGraph.Parent.Parent,'CurrentAxes',subGraph);
        else
            subGraph = graph.findobj('Tag', plotVars{i});
            set(graph,'CurrentAxes',subGraph);
        end
        
        % Set latex as default interpreter
        previousInterpreter  = pub_GraphSetInterpreter('latex');
        
        % Check whether variable is present
        hold on;
        if ~ismember(plotVars{i}, plotVars(isKey(irfLevelChg,plotVars))) 
            % It is not a member
            plotMatrix = NaN(length(periods), 1);
            p = plot(periods, plotMatrix(periods,:));
        elseif ismember(plotVars{i}, {'pii', 'Pii', 'piH', 'piF', 'piStar', 'nr', 'nrStar', 'nrPolicy', 'nrFree', 'polDef', 'r', 'rRN', 'rFree', 'rNatural', 'rPolicy', 'rGap', 'rGood', 'rBad', 'rGoodOnly', 'mShock'}) ...
                || ~isempty(regexp(plotVars{i}, 'rNb\d\d', 'ONCE'))
            plotMatrix = irfLevelChg( plotVars{i} ) * 4 * 100;
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = 'ann. bps';
        elseif ~isempty(regexp(plotVars{i}, 'psi\w\w\d\d', 'ONCE'))
            plotMatrix = irfLevelChg( plotVars{i} ) * 4 * 100;
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = 'ann. bps';
        elseif ~isempty(regexp(plotVars{i}, 'Term_', 'ONCE'))
            plotMatrix = irfLevelChg( plotVars{i} );
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = 'p.p.';
        elseif ismember(plotVars{i}, {'pac', 'pacDif'})
            plotMatrix = irfLevel( plotVars{i} );
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = 'level';
        elseif ismember(plotVars{i}, {'y_gap', 'fxNominalVar', 'pp', 'def', 'tau', 'tauExp', 'bY', 'gY', 'bbeta'})
            plotMatrix = irfLevelChg( plotVars{i} );
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = 'p.p.';
        elseif sum(isinf(irfPercChg( plotVars{i} ))) > 0 %check if it is infinite (ss = 0)
            plotMatrix = irfLevelChg( plotVars{i} );
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = 'level';
        else
            plotMatrix = irfPercChg( plotVars{i} );
            p = plot(periods, plotMatrix(periods,:));
            p(1).Parent.YLabel.String = '\%';
        end
        if setFirstPeriodToNaN
            for iLine = 1:length(p)
                p(iLine).YData(1) = nan;
            end
        end
        
        % Customize lines
        if isempty(customizeLines)
            customizeLines = struct();
            customizeLines.custom_LineWidth = {1, 1, 1, 1, 1, 1, 1};
            customizeLines.custom_LineStyle = {'-', ':', '--', '-.', '-', ':', '--'};
            customizeLines.custom_LineColor = {'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red'};
            customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};
        end
        if isstruct(customizeLines)
            for iLine = 1:length(p)
                p(iLine).LineWidth  = customizeLines.custom_LineWidth{iLine};
                p(iLine).LineStyle  = customizeLines.custom_LineStyle{iLine};
                p(iLine).Color      = customizeLines.custom_LineColor{iLine};
                p(iLine).Marker     = customizeLines.custom_LineMarker{iLine};
            end
        end
        
        %plotLineWidth = 0.5 + 0.1 * (length(findobj(gca,'Type','line'))-1);
        %p.LineWidth = plotLineWidth;
        hold off;

        title(plotVarNames{i});
%         if i==1
%             if ~isempty(legendStr)
%                 legend(legendStr, 'Location', 'northoutside');
%             end
%         end
        ax = gca;
        ax.XAxisLocation = 'bottom';
        %ax.XAxisLocation = 'origin';

        % draw zero axis line
        pub_GraphDrawZeroAxis(p);

        % Adjust y axis
        graphGapPct = 0.1;
        yData = horzcat(p(1).Parent.Children.YData); % concatenate all series
        if min([yData, 0]) == max([yData, 0])
            ax.YLim = [-1 1];
        elseif max(abs(yData)) < tol
            ax.YLim = [-1 1];
        elseif sum(isfinite(max(abs(yData)))) == 0
            ax.YLim = [-1 1];
        else
            ax.YLim = [min([yData, 0]) max([yData, 0])] * (1+graphGapPct);
        end

        if ~isempty(legendStr)
            if legendInSubPlot == 0
                legend(legendStr);
            elseif i == 2
                legendAxes = p;
            end
        end
    end
    
    if legendInSubPlot > 0
        pub_GraphPutLegendInFirstSubplot(...
            subGraph_legend, legendAxes, legendStr, ...
            'Interpreter', 'latex', ...
            'Background', 'boxoff');
    end
    
    % Add title (this command must be run after the all subplots)
    if ~isempty(suptitleStr)
        pub_Suptitle(suptitleStr);
    end
    
    % Restore interpreter
    pub_GraphSetInterpreter(previousInterpreter);
    
end