% Read lattice parameter from cif file
% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
function [a, b, c] = ilm_read_lat_parm_cif(path)
    str_file = fileread(path);
    str_file = strtrim(strsplit(str_file, '\n'));

    a = fn_extract_pattern(str_file, '_cell_length_a', '%s%f%s', 2);
    b = fn_extract_pattern(str_file, '_cell_length_b', '%s%f%s', 2);
    c = fn_extract_pattern(str_file, '_cell_length_c', '%s%f%s', 2);
end

function[str_out] =fn_extract_pattern(str_file, str_ss, patt, idx)
    if(nargin<3)
        idx = 0;
    end
    
    str_iy = str_file{contains(str_file, str_ss)};
    ix = strfind(str_iy, str_ss);
    str_iy = str_iy(ix:end);
    C = textscan(str_iy, patt);
    if(idx==0)
        str_out = C;
    else
        str_out = C{idx};
    end
end