function [atoms, lx, ly, lz, a, b, c, dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    % https://www.webelements.com/gold/crystal_structure.html
    Z = ilm_Z('Au');
    a = 4.0782;
    
    [atoms, lx, ly, lz, a, b, c, dz] = fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end