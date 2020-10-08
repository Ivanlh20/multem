function [atoms, lx, ly, lz, a, b, c, dz] = Primitive_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 4.0780; 
    b = 4.0780; 
    c = 4.0780;
    a = 6.50; 
    b = 6.50;
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 1;
    occ = 1;
    region = 0;
    charge = 0;
    % Au = 79
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [79, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge];
    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end