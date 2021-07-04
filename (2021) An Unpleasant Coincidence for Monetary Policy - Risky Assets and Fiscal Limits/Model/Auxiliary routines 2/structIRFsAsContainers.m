%% Auxiliary functions

function [irfLevel, irfLevelChg, irfPercChg] = structIRFsAsContainers(myirfs, endoNames, ssValues, regimeId)

irfLevel    = containers.Map();
irfLevelChg = containers.Map();
irfPercChg  = containers.Map();

varNames = fieldnames(myirfs);
nVars = length(varNames);
nMdl = size(endoNames, 2);

for ii = 1:nVars
    
    if isa(myirfs.(varNames{ii}), 'ts')
        values = myirfs.(varNames{ii}).values;
    else
        values = myirfs.(varNames{ii}).(['regime_' num2str(regimeId)]).values;
    end

    iSS = cellfun(@(x) find(strcmp(x, varNames{ii})), endoNames);
    irfLevel_aux    = [];
    irfLevelChg_aux = [];
    irfPercChg_aux  = [];
    for iMdl = 1:nMdl
        ss                = ssValues(iSS(iMdl),iMdl);
        irfLevel_aux      = [irfLevel_aux, ss + values(:,iMdl)] ;
        irfLevelChg_aux   = [irfLevelChg_aux, values(:,iMdl)] ;
        irfPercChg_aux    = [irfPercChg_aux, ( (ss + values(:,iMdl)) ./ ss - 1) .* 100] ;
    end
    
    irfLevel(varNames{ii})      = irfLevel_aux;
    irfLevelChg(varNames{ii})   = irfLevelChg_aux;
    irfPercChg(varNames{ii})    = irfPercChg_aux;
    
end

end