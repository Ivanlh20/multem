% Read atomic positions from xyz file
function [atoms, lx, ly, lz] = ilm_read_ap_xyz(path, rmsd_3d)
    ces = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', ...
    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn', 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr', ...
    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',...
    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',...
    'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf'};

    fid = fopen(path,'r');
    fgetl(fid);
    fgetl(fid);
    atoms = []; 
    while feof(fid) == 0
        text = strtrim(upper(sscanf(fgetl(fid),'%c')));
        if (~isempty(text))
            C = textscan(text, '%s %f %f %f');
            atoms = [atoms; [find(strcmpi(ces, strtrim(C{1}))), C{2}, C{3}, C{4}, rmsd_3d, 1.0]];
        end
    end
    fclose(fid);

    atoms(:,2:4) = bsxfun(@minus, atoms(:,2:4), min(atoms(:,2:4)));
    lx = max(atoms(:,2));
    ly = max(atoms(:,3));
    lz = max(atoms(:,4));
end