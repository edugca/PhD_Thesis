function previousInterpreter = pub_GraphSetInterpreter(interpreter)

    previousInterpreter = struct();
    previousInterpreter.defaultAxesTickLabelInterpreter = get(groot, 'defaultAxesTickLabelInterpreter');
    previousInterpreter.defaultLegendInterpreter = get(groot, 'defaultLegendInterpreter');
    previousInterpreter.defaultTextInterpreter = get(groot, 'defaultTextInterpreter');
  
    if isstruct(interpreter)
        set(groot, 'defaultAxesTickLabelInterpreter', interpreter.defaultAxesTickLabelInterpreter);
        set(groot, 'defaultLegendInterpreter', interpreter.defaultLegendInterpreter);
        set(groot, 'defaultTextInterpreter', interpreter.defaultTextInterpreter);
    else
        set(groot, 'defaultAxesTickLabelInterpreter', interpreter);
        set(groot, 'defaultLegendInterpreter', interpreter);
        set(groot, 'defaultTextInterpreter', interpreter);
    end
    
end