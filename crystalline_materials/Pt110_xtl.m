function [atoms, lx, ly, lz, a, b, c, dz] = Pt110_xtl(na, nb, nc, ncu, rms3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    a = 3.9242/sqrt(2); 
    b = 3.9242; 
    c = 3.9242/sqrt(2);
    xtl_parm.a = a;
    xtl_parm.b = b;
    xtl_parm.c = c;
    xtl_parm.nuLayer = 2;
    occ = 1;
    region = 0;
    charge = 0;
    % Pt = 78
    % Z charge x y z rms3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [78, 0.00, 0.00, 0.00, rms3d, 1.0, charge];
    xtl_parm.uLayer(2).atoms = [78, 0.50, 0.50, 0.50, rms3d, 1.0, charge];
    atoms = il_crystal_by_lays(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a; ly = nb*xtl_parm.b; lz = nc*xtl_parm.c;
end