% Plot median, min, max
if config_plotFiltersRangeAllCountries
    
    % Select figure
    if strcmp(currentStatus, 'Advanced')
        f = fAdvanced;
        figure(fAdvanced);
        tl = tlAdvanced;
        legend_axes = legend_axes_advanced;
    elseif strcmp(currentStatus, 'Emerging')
        f = fEmerging;
        figure(fEmerging);
        tl = tlEmerging;
        legend_axes = legend_axes_emerging;
    end
    tile = nexttile();
    
    previousInterpreter = pub_GraphSetInterpreter('latex');
    for iSpec = 1:nSpecs
        % Get specification results
        ts_Filters = datasetFilters.(['Spec_' num2str(iSpec)]);

        %subplot(nLins, nCols, iCountry);
        strTitle = specNames{iSpec};
        p = plot( ...
            ts_Filters.Properties.RowTimes, ...
            ts_Filters.trend_Median, ...
            'black', ...
            'DisplayName', 'Filters median');
        dateaxis('x', 11); % yy
        hold on;
        xArea = [ts_Filters.Properties.RowTimes; flipud(ts_Filters.Properties.RowTimes)];
        yArea = [ts_Filters.trend_Min; flipud(ts_Filters.trend_Max)];
        p2 = fill(xArea, yArea, 'b', 'FaceAlpha', 0.3, 'DisplayName', 'Filters range');
        hold off;
        set(gca, 'FontSize', 14);
        title(strTitle, 'Interpreter', 'latex');
        
        if sum(tl.GridSize) == 3 && legendInFirstSubPlot
            unqLegHands = legendUnq([p p2]);
            strLegend = get([p p2], 'DisplayName')';
            pub_GraphPutLegendInFirstSubplot(...
                legend_axes, [p p2], strLegend, ...
                'Interpreter', 'latex', ...
                'Background', 'boxoff', ...
                'Location', 'northwest');
        end
    end
    pub_GraphSetInterpreter(previousInterpreter);
end