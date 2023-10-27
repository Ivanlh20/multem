% Read lattice parameter from cif file
function [alpha, beta, gamma] = ilm_read_angles_cif(path)
    str_file = fileread(path);
    str_file = strtrim(strsplit(str_file, '\n'));

    alpha = fn_extract_pattern(str_file, '_cell_angle_alpha', '%s%f%s', 2);
    beta = fn_extract_pattern(str_file, '_cell_angle_beta', '%s%f%s', 2);
    gamma = fn_extract_pattern(str_file, '_cell_angle_gamma', '%s%f%s', 2);
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