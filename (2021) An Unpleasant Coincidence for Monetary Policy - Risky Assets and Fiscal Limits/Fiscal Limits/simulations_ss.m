%% Simulation 1: Good TFP

%if ~ismember('epsA', ssExcluded)

    disp('Simulating steady-state values for TFP');

    changeType = 'ss'; % ss or params
    changeName = 'aTildeSS';
    changeLabel = '$\overline{A}$';
    changeLabelLowValue = '-1\%';
    changeLabelHighValue = '+1\%';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName)*0.99;
    paramHigh   = simStruct.(changeType).(changeName)*1.01;

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

%end

%% Simulation 2: Bad TFP

if ~ismember('epsR', shocksExcluded)

    disp('Simulating steady-state values for OmegaTilde');

    changeType = 'ss'; % ss or params
    changeName = 'recTildeSS';
    changeLabel = '$\overline{\tilde{\omega}}$';
    changeLabelLowValue = '-5\%';
    changeLabelHighValue = '+5\%';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName)*0.95;
    paramHigh   = simStruct.(changeType).(changeName)*1.05;

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
    
%% Simulation 3: Default risk

if ~ismember('epsD', shocksExcluded)

    disp('Simulating steady-state values for Def');

    changeType = 'ss'; % ss or params
    changeName = 'defSS';
    changeLabel = '$\overline{\mathcal{D}^r}$';
    changeLabelLowValue = '-5 p.p.';
    changeLabelHighValue = '+5 p.p.';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.05;
    paramHigh   = simStruct.(changeType).(changeName) + 0.05;

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
    
%% Simulation 4: Government expenses

%if ~ismember('epsG', shocksExcluded)

    disp('Simulating steady-state values for G');

    changeType = 'ss'; % ss or params
    changeName = 'gSS';
    changeLabel = '$\overline{G}$';
    changeLabelLowValue = '-1\% $\overline{Y}$';
    changeLabelHighValue = '+1\% $\overline{Y}$';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = (simStruct.(changeType).(changeName) / simStruct.ss.ySS - 0.01) * simStruct.ss.ySS;
    paramHigh   = (simStruct.(changeType).(changeName) / simStruct.ss.ySS + 0.01) * simStruct.ss.ySS;

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
    
%end

%% Simulation 5: Beta

%if ~ismember('epsBeta', shocksExcluded)

    disp('Simulating steady-state values for $\beta$');

    changeType = 'params'; % ss or params
    changeName = 'beta';
    changeLabel = '$\overline{\beta}$';
    changeLabelLowValue = '-0.25 p.p.';
    changeLabelHighValue = '+0.25 p.p.';
    simGraphName = ['fiscalLimitDist_', changeName, '.png'];
    simGraphName_CDF = ['fiscalLimitCDF_', changeName, '.png'];

    paramLow    = simStruct.(changeType).(changeName) - 0.0025;
    paramHigh   = simStruct.(changeType).(changeName) + 0.0025;

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

%end

%% Consolidate plots

simGraphName_Consol = ['fiscalLimitPDFandCDFAll_SS.png'];

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
