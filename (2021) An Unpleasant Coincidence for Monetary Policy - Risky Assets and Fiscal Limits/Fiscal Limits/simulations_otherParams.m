%% Simulation 1: gammaGPSI

if ~ismember('epsG', shocksExcluded)

    disp('Simulating steady-state values for TFP');

    changeType = 'params'; % ss or params
    changeName = 'gammaGPSI';
    changeLabel = '$\gamma_{G\Psi}$';
    changeLabelLowValue = '0.0';
    changeLabelHighValue = '0.3';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = 0.0;
    paramHigh   = 0.3;

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

%% Simulation 2: alphaG

if ~ismember('epsG', shocksExcluded)

    disp('Simulating steady-state values for TFP');

    changeType = 'params'; % ss or params
    changeName = 'alphaG';
    changeLabel = '$\alpha_G$';
    changeLabelLowValue = '-1.0';
    changeLabelHighValue = '1.0';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = -1 ;
    paramHigh   = 1 ;

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

%% Simulation 3: fracNR

if ~ismember('epsG', shocksExcluded)

    disp('Simulating steady-state values for TFP');

    changeType = 'params'; % ss or params
    changeName = 'fracNR';
    changeLabel = '$\gamma^{NR}$';
    changeLabelLowValue = '0.0';
    changeLabelHighValue = '0.8';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    %paramLow    = simStruct.(changeType).(changeName) - 0.4;
    %paramHigh   = simStruct.(changeType).(changeName) + 0.5;
    paramLow    = 0.0;
    paramHigh   = 0.8;

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

simGraphName_Consol = ['fiscalLimitPDFandCDFAll_otherParams.png'];

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
