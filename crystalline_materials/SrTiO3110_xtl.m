function [atoms, lx, ly, lz, a, b, c, dz] = SrTiO3110_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 3.9050; 
    b = 3.9050*sqrt(2); 
    c = 3.9050*sqrt(2);
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 4;
    occ = 1;
    region = 0;
    charge = 0;
    % Sr = 38, Ti = 22; O = 8
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [38, 0.50, 0.00, 0.00, rmsd_3d, occ, region, charge; 8, 0.50, 0.50, 0.00, rmsd_3d, occ, region, charge; 22, 0.0, 0.5, 0.0, rmsd_3d, occ, region, charge];
    xtl_parm.uLayer(2).atoms = [8, 0.00, 0.25, 0.25, rmsd_3d, occ, region, charge; 8, 0.00, 0.75, 0.25, rmsd_3d, occ, region, charge];
    xtl_parm.uLayer(3).atoms = [8, 0.50, 0.00, 0.50, rmsd_3d, occ, region, charge; 38, 0.50, 0.50, 0.50, rmsd_3d, occ, region, charge; 22, 0.00, 0.00, 0.50, rmsd_3d, occ, region, charge];
    xtl_parm.uLayer(4).atoms = [8, 0.00, 0.25, 0.75, rmsd_3d, occ, region, charge; 8, 0.00, 0.75, 0.75, rmsd_3d, occ, region, charge];

    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end