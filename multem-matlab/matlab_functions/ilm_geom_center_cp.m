function[p0]=ilm_geom_center_cp(atoms)
    p0 = mean(atoms);     % rotation point
    [~, ii] = sort(sqrt(sum((atoms-p0).^2, 2)));
    p0 = atoms(ii(1), :); 
end