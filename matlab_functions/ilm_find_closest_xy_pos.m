function[x_c, y_c]=ilm_find_closest_xy_pos(atoms, x_c, y_c)
    [~, ii] = sort((atoms(:, 2)-x_c).^2+(atoms(:, 3)-y_c).^2);
    x_c = atoms(ii(1), 2);
    y_c = atoms(ii(1), 3);
end