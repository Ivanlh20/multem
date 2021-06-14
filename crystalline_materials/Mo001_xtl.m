function [atoms, lx, ly, lz, a, b, c, dz] = Mo001_xtl(na, nb, nc, ncu, rmsd_3d)
    if(nargin<5)
        rmsd_3d = 0.085;
    end
    
    Z = ilm_Z('Mo');
    a = 3.147;
    
    [atoms, lx, ly, lz, a, b, c, dz] = bcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d);
end