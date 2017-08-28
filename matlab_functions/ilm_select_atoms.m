function[atoms]=ilm_select_atoms(atoms, x_min, x_max, y_min, y_max)
    x = atoms(:, 2);
    y = atoms(:, 3);    
    atoms = atoms((x_min<=x)&(x<x_max)&(y_min<=y)&(y<y_max), :);   
end