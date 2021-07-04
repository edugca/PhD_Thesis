%% Simulation 5: TFP Std

if ~ismember('epsA', shocksExcluded)

    disp('Simulating parameter values for TFP Std');

    changeType = 'params'; % ss or params
    changeName = 'sigmaA';
    changeLabel = '$\sigma_A$';
    changeLabelLowValue = '0.1 $\overline{\sigma_A}$';
    changeLabelHighValue = '2.0 $\overline{\sigma_A}$';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName)*0.1;
    paramHigh   = simStruct.(changeType).(changeName)*2.0;

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

%% Simulation 6: recTilde Std

if ~ismember('epsR', shocksExcluded)

    disp('Simulating parameter values for recTilde Std');

    changeType = 'params'; % ss or params
    changeName = 'sigmaR';
    changeLabel = '$\sigma_\omega$';
    changeLabelLowValue = '0.1 $\overline{\sigma_\omega}$';
    changeLabelHighValue = '2.0 $\overline{\sigma_\omega}$';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName)*0.1;
    paramHigh   = simStruct.(changeType).(changeName)*2.0;

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
    
%% Simulation 7: def Std

if ~ismember('epsD', shocksExcluded)

    disp('Simulating parameter values for def Std');

    changeType = 'params'; % ss or params
    changeName = 'sigmaD';
    changeLabel = '$\sigma_D$';
    changeLabelLowValue = '0.1 $\overline{\sigma_D}$';
    changeLabelHighValue = '2.0 $\overline{\sigma_D}$';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName)*0.1;
    paramHigh   = simStruct.(changeType).(changeName)*2.0;

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

%% Simulation 8: G Std

if ~ismember('epsG', shocksExcluded)

    disp('Simulating parameter values for G Std');

    changeType = 'params'; % ss or params
    changeName = 'sigmaG';
    changeLabel = '$\sigma_G$';
    changeLabelLowValue = '0.1 $\overline{\sigma_G}$';
    changeLabelHighValue = '2.0 $\overline{\sigma_G}$';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName)*0.1;
    paramHigh   = simStruct.(changeType).(changeName)*2.0;

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

%% Simulation 9: Beta Std

if ~ismember('epsBeta', shocksExcluded)

    disp('Simulating parameter values for Beta Std');

    changeType = 'params'; % ss or params
    changeName = 'sigmaBeta';
    changeLabel = '$\sigma_\beta$';
    changeLabelLowValue = '0';
    changeLabelHighValue = '$\overline{\sigma_\beta}$ + 0.001';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = 0; % it can't be negative
    paramHigh   = simStruct.(changeType).(changeName) + 0.001;

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

simGraphName_Consol = ['fiscalLimitPDFandCDFAll_Sigma.png'];

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