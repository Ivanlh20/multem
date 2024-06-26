function [atoms, lx, ly, lz, a, b, c, dz] = Cu001_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    % https://www.webelements.com/copper/crystal_structure.html
    Z = ilm_Z('Cu');
    a = 3.6149;
    
    [atoms, lx, ly, lz, a, b, c, dz] = fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end