function [atoms] = center_specimen(atoms, lx, ly, lz)
atoms(:, 2:4) = bsxfun(@minus, atoms(:, 2:4), min(atoms(:, 2:4)));
lxs = 0.5*(lx-max(atoms(:, 2)));
lys = 0.5*(ly-max(atoms(:, 3)));
lzs = 0.5*(lz-max(atoms(:, 4)));
atoms(:, 2:4) = bsxfun(@plus, atoms(:, 2:4), [lxs, lys, lzs]);