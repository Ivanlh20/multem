function[p0] = ilm_geom_center_cp(atoms)
    p0 = ilm_find_closest_atom_pos(mean(atoms));
end