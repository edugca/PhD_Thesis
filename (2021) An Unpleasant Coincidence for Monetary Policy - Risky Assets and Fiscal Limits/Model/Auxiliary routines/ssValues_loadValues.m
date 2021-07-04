%% Steady state values: load values

ssValues = fieldnames(simStruct.ss);

for ii=1:length(ssValues)
   eval([ssValues{ii}, ' = simStruct.ss.', ssValues{ii}, ';']); 
end