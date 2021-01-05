% Read atomic positions from cif file
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
function [atoms, lx, ly, lz] = ilm_read_ap_cif(path, rmsd_3d_0, pbc, na, nb, nc)
    if(nargin<6)
        nc = 1;
    end
    
    if(nargin<5)
        nb = 1;
    end
    
     if(nargin<4)
        na = 1;
    end

    % for anisotropic atomic displacement check out the following website
    % https://www.iucr.org/__data/iucr/cif/software/oostar/oostar_ddl1/bin/data/cifdic.c91.list.html
    if(nargin<3)
        pbc = true;
    end
    
    if(nargin<2)
        rmsd_3d_0 = 0.085;
    end
    
    occ_0 = 1;

    str_file = fileread(path);
    str_file = strtrim(strsplit(str_file, '\n'));

    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;

    xtl_parm.a = fn_extract_pattern(str_file, '_cell_length_a', '%s%f%s', 2);
    xtl_parm.b = fn_extract_pattern(str_file, '_cell_length_b', '%s%f%s', 2);
    xtl_parm.c = fn_extract_pattern(str_file, '_cell_length_c', '%s%f%s', 2);

    xtl_parm.alpha = fn_extract_pattern(str_file, '_cell_angle_alpha', '%s%f%s', 2);
    xtl_parm.beta = fn_extract_pattern(str_file, '_cell_angle_beta', '%s%f%s', 2);
    xtl_parm.gamma = fn_extract_pattern(str_file, '_cell_angle_gamma', '%s%f%s', 2);

    % https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_IT_number.html
    if(fn_get_idx(str_file, '_space_group_IT_number')>0)
        xtl_parm.sgn = fn_extract_pattern(str_file, '_space_group_IT_number', '%s%d%s', 2);
    else
        xtl_parm.sgn = fn_extract_pattern(str_file, '_symmetry_Int_Tables_number', '%s%d%s', 2);
    end
    
    xtl_parm.pbc = pbc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iy_0 = find(contains(str_file, '_atom_site_fract_x'), 1, 'first');
    iy_0 = iy_0-1;
    while ~contains(str_file{iy_0}, 'loop_')
        iy_0 = iy_0-1;
    end
    iy_0 = iy_0 + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iy_e = iy_0+1;
    while contains(str_file{iy_e}, '_atom_site_')
        iy_e = iy_e+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_col = (iy_e-iy_0);
    map = containers.Map(iy_0:(iy_e-1), 1:n_col);
    
    iy_Z = fn_get_idx(str_file, '_atom_site_type_symbol');
    iy_x = fn_get_idx(str_file, '_atom_site_fract_x');
    iy_y = fn_get_idx(str_file, '_atom_site_fract_y');
    iy_z = fn_get_idx(str_file, '_atom_site_fract_z');
    
    iy_rmsd_3d = fn_get_idx(str_file, '_atom_site_B_iso_or_equiv');
    iy_occ = fn_get_idx(str_file, '_atom_site_occupancy');
    str_patt = repmat('%s', 1, n_col);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    asym_uc = [];
    n_rows = length(str_file);
    iy = iy_e;
    while true
        if((iy>n_rows) || strcmp(str_file{iy}, ''))
            break;
        end
        
        str_text = textscan(str_file{iy}, str_patt);
        if(fn_break(str_text{map(iy_Z)}))
            break;
        end
        
        Z = fn_read_Z(str_text{map(iy_Z)}{1});
        x = fn_read_number(str_text{map(iy_x)}{1});
        y = fn_read_number(str_text{map(iy_y)}{1});
        z = fn_read_number(str_text{map(iy_z)}{1});
        
         % https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_B_iso_or_equiv.html
         if(ilm_chk_bound(iy_rmsd_3d, iy_0, iy_e))
            rmsd_3d = fn_read_rmsd_3d(str_text{map(iy_rmsd_3d)}{1}, rmsd_3d_0);
         else
            rmsd_3d = rmsd_3d_0;
         end
         
         % https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_occupancy.html
         if(ilm_chk_bound(iy_occ, iy_0, iy_e))
            occ = fn_read_rmsd_3d(str_text{map(iy_occ)}{1}, occ_0);
         else
            occ = occ_0;
         end
         
        asym_uc = [asym_uc;Z, x, y, z, rmsd_3d, occ];
        iy = iy + 1;
    end

    xtl_parm.asym_uc = asym_uc;
    xtl_parm.base = [];
    
    atoms = ilc_xtl_build(xtl_parm);
    
    xyz = atoms(:, 2:4);
    xyz = xyz - min(xyz);
    [lx, ly, lz] = ilm_vect_assign(max(xyz));
    atoms(:, 2:4) = xyz;
    
    [lx, ly, lz] = ilm_vect_assign([max(lx, xtl_parm.a), max(ly, xtl_parm.b), max(lz, xtl_parm.c)]);
end

function[bb] =fn_break(str_ss)
    if(size(str_ss, 1)>0)
        str_ss = str_ss{1};
        bb = ilm_Z(str_ss(isletter(str_ss)))<1;
    else
        bb = true;
    end
end

function[Z] =fn_read_Z(str_ss)
    Z = ilm_Z(str_ss(isletter(str_ss)));
end

function[n] =fn_read_number(str_ss, n_0)
    if(nargin<2)
        n_0 = 0;
    end
    
    str_ss = eraseBetween(str_ss, '(', ')', 'Boundaries', 'inclusive');
    n = str2double(strtrim(str_ss));
    
    if(isnan(n))
        n = n_0;
    end
end

function[rmsd_3d] =fn_read_rmsd_3d(str_ss, rmsd_3d_0)
    B = fn_read_number(str_ss);
    if(B==0)
        rmsd_3d = rmsd_3d_0;
    else
        rmsd_3d = sqrt(B/(8*pi^2));
    end
end

function[idx] =fn_get_idx(str_file, str_ss)
    idx = find(contains(str_file, str_ss), 1, 'first');
    if(numel(idx)==0)
        idx = -1;
    end
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