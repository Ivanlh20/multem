function [atoms, lx, ly, lz, a, b, c, dz] = Ni110_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    % https://www.webelements.com/nickel/crystal_structure.html
    Z = ilm_Z('Ni');
    a = 3.5240;
    
    [atoms, lx, ly, lz, a, b, c, dz] = fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end