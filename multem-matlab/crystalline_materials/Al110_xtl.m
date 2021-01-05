function [atoms, lx, ly, lz, a, b, c, dz] = Al110_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    % https://www.webelements.com/aluminium/crystal_structure.html
    Z = ilm_Z('Al');
    a = 4.0495;
    
    [atoms, lx, ly, lz, a, b, c, dz] = fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end