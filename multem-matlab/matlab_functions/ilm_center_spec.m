function [atoms] = ilm_center_spec(atoms, lx, ly, lz)
    atoms(:, 2:4) = bsxfun(@minus, atoms(:, 2:4), min(atoms(:, 2:4)));
    atoms(:, 2:4) = bsxfun(@plus, atoms(:, 2:4), 0.5*([lx, ly, lz]-max(atoms(:, 2:4))));
end