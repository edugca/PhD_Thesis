%% Draw axis at y = 0
function pub_GraphDrawZeroAxis(graph)

if ismember('XData', fieldnames(graph))
    xValues = graph.XData;
    
    originLine = line(...
        [xValues(1) xValues(end)], ...
        [0 0], ...
        'Color','black','LineStyle','--','HandleVisibility','off');
else
    xValues = graph.XLim;
    
    originLine = line(...
        graph, ...
        [xValues(1) xValues(end)], ...
        [0 0], ...
        'Color','black','LineStyle','--','HandleVisibility','off');
end




end