function r = getBBOBRootPath(rootName)
    r = pwd;
    ind = max(strfind(r,[filesep rootName filesep]));
    r(ind+length(rootName)+1:end) = [];
end