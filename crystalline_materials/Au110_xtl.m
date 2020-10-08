function [atoms, lx, ly, lz, a, b, c, dz] = Au110_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 4.0780/sqrt(2); 
    b = 4.0780; 
    c = 4.0780/sqrt(2);
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 2;
    occ = 1;
    region = 0;
    charge = 0;
    % Au = 79
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [79, 0.00, 0.00, 0.00, rmsd_3d, 1.0];
    xtl_parm.uLayer(2).atoms = [79, 0.50, 0.50, 0.50, rmsd_3d, 1.0];
    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end