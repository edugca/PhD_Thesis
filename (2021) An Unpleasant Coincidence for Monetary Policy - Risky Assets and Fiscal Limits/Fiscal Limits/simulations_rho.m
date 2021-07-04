%% Simulation 5: TFP AR

if ~ismember('epsA', shocksExcluded)

    disp('Simulating parameter values for TFP AR');

    changeType = 'params'; % ss or params
    changeName = 'rhoA';
    changeLabel = '$\rho_A$';
    changeLabelLowValue     = '-0.025';
    changeLabelHighValue    = '+0.025';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.025;

    plotFiscalLimits;

    iTFL = iTFL + 1;
    tFL{iTFL, 1} = simStruct.(changeType).(changeName);
    tFL{iTFL, 2} = paramLow;
    tFL{iTFL, 3} = paramHigh;
    tFL{iTFL, 4} = median(fiscalLimit_1);
    tFL{iTFL, 5} = median(fiscalLimit_2);
    tFL{iTFL, 6} = median(fiscalLimit_3);
    tFL{iTFL, 7} = std(fiscalLimit_1);
    tFL{iTFL, 8} = std(fiscalLimit_2);
    tFL{iTFL, 9} = std(fiscalLimit_3);

end

%% Simulation 6: recTilde AR

if ~ismember('epsR', shocksExcluded)

    disp('Simulating parameter values for recTilde AR');

    changeType = 'params'; % ss or params
    changeName = 'rhoR';
    changeLabel = '$\rho_\omega$';
    changeLabelLowValue     = '-0.025';
    changeLabelHighValue    = '+0.025';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.025;

    plotFiscalLimits;

    iTFL = iTFL + 1;
    tFL{iTFL, 1} = simStruct.(changeType).(changeName);
    tFL{iTFL, 2} = paramLow;
    tFL{iTFL, 3} = paramHigh;
    tFL{iTFL, 4} = median(fiscalLimit_1);
    tFL{iTFL, 5} = median(fiscalLimit_2);
    tFL{iTFL, 6} = median(fiscalLimit_3);
    tFL{iTFL, 7} = std(fiscalLimit_1);
    tFL{iTFL, 8} = std(fiscalLimit_2);
    tFL{iTFL, 9} = std(fiscalLimit_3);
    
end

%% Simulation 7: def AR

if ~ismember('epsD', shocksExcluded)

    disp('Simulating parameter values for def AR');

    changeType = 'params'; % ss or params
    changeName = 'rhoD';
    changeLabel = '$\rho_D$';
    changeLabelLowValue     = '-0.025';
    changeLabelHighValue    = '+0.025';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.025;

    plotFiscalLimits;

    iTFL = iTFL + 1;
    tFL{iTFL, 1} = simStruct.(changeType).(changeName);
    tFL{iTFL, 2} = paramLow;
    tFL{iTFL, 3} = paramHigh;
    tFL{iTFL, 4} = median(fiscalLimit_1);
    tFL{iTFL, 5} = median(fiscalLimit_2);
    tFL{iTFL, 6} = median(fiscalLimit_3);
    tFL{iTFL, 7} = std(fiscalLimit_1);
    tFL{iTFL, 8} = std(fiscalLimit_2);
    tFL{iTFL, 9} = std(fiscalLimit_3);

end
    
%% Simulation 8: G AR

if ~ismember('epsG', shocksExcluded)

    disp('Simulating parameter values for GG AR');

    changeType = 'params'; % ss or params
    changeName = 'rhoGG';
    changeLabel = '$\rho_{GG}$';
    changeLabelLowValue     = '-0.025';
    changeLabelHighValue    = '+0.025';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.025;

    plotFiscalLimits;

    iTFL = iTFL + 1;
    tFL{iTFL, 1} = simStruct.(changeType).(changeName);
    tFL{iTFL, 2} = paramLow;
    tFL{iTFL, 3} = paramHigh;
    tFL{iTFL, 4} = median(fiscalLimit_1);
    tFL{iTFL, 5} = median(fiscalLimit_2);
    tFL{iTFL, 6} = median(fiscalLimit_3);
    tFL{iTFL, 7} = std(fiscalLimit_1);
    tFL{iTFL, 8} = std(fiscalLimit_2);
    tFL{iTFL, 9} = std(fiscalLimit_3);

end
    
%% Simulation 8: G AR

if ~ismember('epsG', shocksExcluded)

    disp('Simulating parameter values for GY AR');

    changeType = 'params'; % ss or params
    changeName = 'rhoGY';
    changeLabel = '$\rho_{GY}$';
    changeLabelLowValue     = '-0.025';
    changeLabelHighValue    = '+0.025';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.025;

    plotFiscalLimits;

    iTFL = iTFL + 1;
    tFL{iTFL, 1} = simStruct.(changeType).(changeName);
    tFL{iTFL, 2} = paramLow;
    tFL{iTFL, 3} = paramHigh;
    tFL{iTFL, 4} = median(fiscalLimit_1);
    tFL{iTFL, 5} = median(fiscalLimit_2);
    tFL{iTFL, 6} = median(fiscalLimit_3);
    tFL{iTFL, 7} = std(fiscalLimit_1);
    tFL{iTFL, 8} = std(fiscalLimit_2);
    tFL{iTFL, 9} = std(fiscalLimit_3);
    
end

%% Simulation 9: Beta AR

if ~ismember('epsBeta', shocksExcluded)

    disp('Simulating parameter values for Beta AR');

    changeType = 'params'; % ss or params
    changeName = 'rhoBeta';
    changeLabel = '$\rho_{\beta}$';
    changeLabelLowValue     = '-0.025';
    changeLabelHighValue    = '+0.025';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.025;

    plotFiscalLimits;

    iTFL = iTFL + 1;
    tFL{iTFL, 1} = simStruct.(changeType).(changeName);
    tFL{iTFL, 2} = paramLow;
    tFL{iTFL, 3} = paramHigh;
    tFL{iTFL, 4} = median(fiscalLimit_1);
    tFL{iTFL, 5} = median(fiscalLimit_2);
    tFL{iTFL, 6} = median(fiscalLimit_3);
    tFL{iTFL, 7} = std(fiscalLimit_1);
    tFL{iTFL, 8} = std(fiscalLimit_2);
    tFL{iTFL, 9} = std(fiscalLimit_3);
    
end

%% Consolidate plots

simGraphName_Consol = ['fiscalLimitPDFandCDFAll_Rho.png'];

figHandles = flipud(get(groot, 'Children'));
[consolidatedPlot, subPlotAxes] = edu_TransformFiguresInSubPlots(figHandles,2);
consolidatedPlot;
edu_GraphSetInterpreter('latex');
for iAx = 1:2:(length(subPlotAxes)-1)
    legend(subPlotAxes(iAx), 'boxoff');
    legend(subPlotAxes(iAx), 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 3);
end
xlabel(subPlotAxes(end-1:end), '$\frac{\overline{\mathcal{B}}}{4\overline{Y}} \, (\%)$');
ylabel(subPlotAxes(1:2:(end-1)), 'probability density', 'Interpreter', 'latex');
ylabel(subPlotAxes(2:2:end), 'cumm. prob. (\%)', 'Interpreter', 'latex');
set(subPlotAxes, 'XTick', 50:20:250);
set(subPlotAxes, 'XGrid', 'on');
set(subPlotAxes, 'YGrid', 'on');
linkaxes(subPlotAxes, 'x');
subPlotAxes(1).XLim = [50, 150];
%linkaxes(subPlotAxes(1:2:(end-1)), 'y');
linkaxes(subPlotAxes(2:2:end), 'y');
%edu_Suptitle('CDFs and PDFs of Fiscal Limits','FontSize', 24);
set(findobj(gcf,'type','axes'),'FontSize', 14);grid on

set(gcf, 'Position',  [100, 100, 800, length(figHandles)/2 * 300]); % resize figure
saveas(consolidatedPlot,[imagesFolder, simGraphName_Consol]);