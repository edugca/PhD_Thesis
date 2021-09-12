function pub_GraphPutLegendInFirstSubplot(firstSubPlotAxis, legendAxes, ...
    legendStrings, varargin)

% default
fontSize    = 16;
interpreter = '';
location    = 'best';
background      = 'boxon';

for ii = 1:2:length(varargin)
    if strcmp(varargin(ii), 'FontSize')
       fontSize = varargin{ii+1};
    elseif strcmp(varargin(ii), 'Interpreter')
       interpreter = varargin{ii+1};
    elseif strcmp(varargin(ii), 'Location')
       location = varargin{ii+1};
    elseif strcmp(varargin(ii), 'Background')
       background = varargin{ii+1};
    end
end

% Add single legend slot
poshL = get(firstSubPlotAxis,'position');     % Getting its position
lgd = legend(firstSubPlotAxis, ...
             legendAxes, ...
             legendStrings, ...
        'FontSize', fontSize, ...
        'Interpreter', interpreter, ...
        'Location', location);
set(lgd,'position', poshL);      % Adjusting legend's position
axis(firstSubPlotAxis, 'off');                 % Turning its axis off
legend(firstSubPlotAxis, background);

end