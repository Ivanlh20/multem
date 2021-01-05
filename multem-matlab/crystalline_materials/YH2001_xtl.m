function [atoms, lx, ly, lz, a, b, c, dz] = YH2001_xtl(na, nb, nc, ncu, rmsd_3d_Y, rmsd_3d_H)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 5.2000; 
    b = 5.2000; 
    c = 5.2000; 
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 4;
    occ = 1;
    region = 0;
    charge = 0;
    % Au = 79
    % Z x y z rmsd_3d occupancy charge
    xtl_parm.uLayer(1).atoms = [39, 0.0, 0.0, 0.0, rmsd_3d_Y, occ, region, charge; 39, 0.5, 0.5, 0.0, rmsd_3d_Y, occ, region, charge];
    xtl_parm.uLayer(2).atoms = [1, 0.25, 0.25, 0.25, rmsd_3d_H, occ, region, charge; 1, 0.75, 0.25, 0.25, rmsd_3d_H, occ, region, charge; 1, 0.25, 0.75, 0.25, rmsd_3d_H, occ, region, charge; 1, 0.75, 0.75, 0.25, rmsd_3d_H, occ, region, charge];
    xtl_parm.uLayer(3).atoms = [39, 0.0, 0.5, 0.5, rmsd_3d_Y, occ, region, charge; 39, 0.5, 0.0, 0.5, rmsd_3d_Y, occ, region, charge];
    xtl_parm.uLayer(4).atoms = [1, 0.25, 0.25, 0.75, rmsd_3d_H, occ, region, charge; 1, 0.75, 0.25, 0.75, rmsd_3d_H, occ, region, charge; 1, 0.25, 0.75, 0.75, rmsd_3d_H, occ, region, charge; 1, 0.75, 0.75, 0.75, rmsd_3d_H, occ, region, charge];

    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end