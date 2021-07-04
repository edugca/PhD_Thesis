%% Parameter values: load values

params = fieldnames(simStruct.params);

for ii=1:length(params)
   eval([params{ii}, ' = simStruct.params.', params{ii}, ';']); 
end