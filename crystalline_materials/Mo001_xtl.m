function [atoms, lx, ly, lz, a, b, c, dz] = Mo001_xtl(na, nb, nc, ncu, rms3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 3.1470; 
    b = 3.1470; 
    c = 3.1470;
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 2;
    occ = 1;
    region = 0;
    charge = 0;
    % Mo = 42
    % Z x y z rms3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [42, 0.0, 0.0, 0.0, rms3d, occ, region, charge];
    xtl_parm.uLayer(2).atoms = [42, 0.5, 0.5, 0.5, rms3d, occ, region, charge];
    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end