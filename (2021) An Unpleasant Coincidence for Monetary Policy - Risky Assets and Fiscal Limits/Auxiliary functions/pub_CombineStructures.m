function S = pub_CombineStructures(dim, varargin)
    % specificProperty = [] ; caoncatenates the values of the structure

    if length(varargin) == 1
        varargin = cell2mat(varargin);
        aux = varargin;
        varargin = {};
        for iStruct = 1:length(aux)
            varargin{iStruct} = aux(iStruct);
        end
    end

    % Check whether there is only one structure
    if length(varargin) == 1
        S = varargin{1};
        return 
    end
    
    F = cellfun(@fieldnames,varargin,'uni',0);
    
    assert(isequal(F{:}),'All structures must have the same field names.')
    
    T = [varargin{:}];
    S = struct();
    F = F{1};
    
    for k = 1:numel(F)
        S.(F{k}) = cat(dim,T.(F{k}));
    end
    
end