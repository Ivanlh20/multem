function [atoms, lx, ly, lz, a, b, c, dz] = Pt110_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    % https://www.webelements.com/platinum/crystal_structure.html
    Z = ilm_Z('Pt');
    a = 3.9242;
    
    [atoms, lx, ly, lz, a, b, c, dz] = fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end