function [atoms, lx, ly, lz, a, b, c, dz] = SrTiO3001_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 3.9050; 
    b = 3.9050; 
    c = 3.9050;
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 2;
    occ = 1;
    region = 0;
    charge = 0;
    % Sr = 38, Ti = 22; O = 8
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [38, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge; 8, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge];
    xtl_parm.uLayer(2).atoms = [8, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge; 8, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge; 22, 0.5, 0.5, 0.5, rmsd_3d, occ, region, charge];
    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end