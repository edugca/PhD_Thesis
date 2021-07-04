%% No government list

if riseConfigVector(refMdl).conf_hasGovernment == false
    
    %%%%%%%%%%%
    varlist         = { 'aTilde', 'recTilde', 'def', 'mShock', ...
                        'y', 'c', 'w', 'n', ...
                        'nrPolicy' , 'rRN', 'rGap', 'Pii'};
    varlist         = varlist(find(ismember(varlist, endo_names{refMdl}))); % clear variables not present
    vlocs           = locate_variables(varlist, endo_names{refMdl});
    vtexNames       = endo_texnames{refMdl}(vlocs);
    
    plotVars        = varlist;
    plotVarNames    = vtexNames;
    %%%%%%%%%%%
    
    nCols               = 4;
    nLins               = ceil(length(plotVars) / nCols);
    nPlotPeriods        = 1:20;
    suptitleStr         = '';
    legendStr           = setUpNames;
    legendFirstSubPlot  = true;
    setFirstPeriodToNaN = true;
end

if riseConfigVector(refMdl).conf_hasGovernment == true
    
    %%%%%%%%%%%
    varlist         = varlist(find(ismember(varlist, endo_names{refMdl}))); % clear variables not present
    vlocs           = locate_variables(varlist, endo_names{refMdl});
    vtexNames       = endo_texnames{refMdl}(vlocs);
    
    plotVars        = varlist;
    plotVarNames    = vtexNames;
    %%%%%%%%%%%
    
    nCols               = 4;
    nLins               = ceil(length(plotVars) / nCols);
    nPlotPeriods        = 1:40;
    suptitleStr         = '';
    legendStr           = setUpNames;
    legendFirstSubPlot  = true;
    setFirstPeriodToNaN = true;
    
end

numColors = 7;
spaceColors = linspace(0,pi*3,1000);
selColors = linspecer(numColors);

customizeLines = struct();
customizeLines.custom_LineWidth = {1, 2, 1, 1, 1, 1, 1};
customizeLines.custom_LineStyle = {'-', ':', '-.', '--', '-', ':', '--'};
customizeLines.custom_LineColor = {selColors(1,:), selColors(2,:), selColors(3,:), 'blue', 'red', 'red', 'red'};
customizeLines.custom_LineMarker = {'none', 'none', 'none', 'none', 'none', 'none', 'none'};