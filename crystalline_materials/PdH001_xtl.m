function [atoms, lx, ly, lz, a, b, c, dz] = PdH001_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 4.025; 
    b = 4.025; 
    c = 4.025;
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 2;
    occ = 1;
    region = 0;
    charge = 0;
    % Pd = 46, H = 1
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [46, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge; 46, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge;...
                               1, 0.5, 0.0, 0.0, rmsd_3d, occ, region, charge; 1, 0.0, 0.5, 0.0, rmsd_3d, occ, region, charge];
    xtl_parm.uLayer(2).atoms = [46, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge; 46, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge;...
                               1, 0.0, 0.0, 0.5, rmsd_3d, occ, region, charge; 1, 0.5, 0.5, 0.5, rmsd_3d, occ, region, charge];
    atoms = ilc_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end