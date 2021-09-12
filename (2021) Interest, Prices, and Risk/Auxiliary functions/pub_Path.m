function newPath = pub_Path(pathString, prefixPC, prefixMac)

pcFileSep   = '\';
macFileSep  = '/';

if ~exist('prefixPC', 'var')
   prefixPC = ''; 
end

if ~exist('prefixMac', 'var')
   prefixMac = ''; 
end

if ispc
    if ~startsWith(pathString, prefixPC)
        pathString = [prefixPC strrep(pathString, macFileSep, filesep)];
    end
elseif ismac
    if ~startsWith(pathString, prefixMac)
        pathString = [prefixMac strrep(pathString, pcFileSep, filesep)];
    end
end

% Remove double path separator
pathString = strrep(pathString, [filesep filesep], filesep);

newPath = pathString;

end