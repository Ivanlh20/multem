% Read atomic positions from xyz file
function [atoms, lx, ly, lz] = ilm_read_ap_xyz(path, rmsd_3d)
    if(nargin<2)
        rmsd_3d = 0.085;
    end
    
    fid = fopen(path, 'r');
    fgetl(fid);
    fgetl(fid);
    atoms = [];
    while feof(fid) == 0
        text = strtrim(upper(sscanf(fgetl(fid), '%c')));
        if (~isempty(text))
            C = textscan(text, '%s %f %f %f');
            atoms = [atoms;[ilm_Z(C{1}), C{2}, C{3}, C{4}, rmsd_3d, 1.0]];
        end
    end
    fclose(fid);

    xyz = atoms(:, 2:4);
    xyz = xyz - min(xyz);
    [lx, ly, lz] = ilm_vect_assign(max(xyz));
    atoms(:, 2:4) = xyz;
end