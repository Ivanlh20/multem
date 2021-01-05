function[str]=ilm_replace_filesep(str, ss_rep)
    if nargin < 2
        ss_rep = filesep;
    end
    
    if(ispc)
        ss = '/';
    else
        ss = '\';
    end
    
    str = strrep(str, ss, ss_rep);
end