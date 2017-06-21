function[atoms]=extract_atomic_pos(path, ncols, iconf_0)
    sextract = strtrim(repmat('%f ', 1, ncols));
    
    fileID = fopen(path,'r');
    for i1 = 1:3
        fgetl(fileID);
    end
    natoms = str2double(deblank(fgetl(fileID)));
    fclose(fileID);

    function[atoms_t]= read_pos(fileID, sextract, natoms)
        for i2 =1:9
            fgetl(fileID);
        end
        t = textscan(fileID, sextract,'Delimiter','\t', 'MultipleDelimsAsOne', natoms);
        atoms_t = [t{1}, t{2}, t{3}];
    end

    atoms = zeros(natoms, 3);
    fileID = fopen(path, 'r');
    ic = 1;
    while (~feof(fileID) && (ic<=iconf_0))
        atoms_conf = read_pos(fileID, sextract, natoms);
        if(ic==iconf_0)
            atoms = atoms_conf;
            break;
        end
        ic = ic + 1;
    end
    fclose(fileID);
end