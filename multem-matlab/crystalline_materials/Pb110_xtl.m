function [atoms, lx, ly, lz, a, b, c, dz] = Pb110_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    % https://www.webelements.com/lead/crystal_structure.html
    Z = ilm_Z('Pb');
    a = 4.9508;
    
    [atoms, lx, ly, lz, a, b, c, dz] = fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end