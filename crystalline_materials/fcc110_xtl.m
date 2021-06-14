function [atoms, lx, ly, lz, a, b, c, dz] = fcc110_xtl(Z, l_c, na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    
    xtl_parm.a = l_c/sqrt(2);
    xtl_parm.b = l_c;
    xtl_parm.c = l_c/sqrt(2);
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;

    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];

    occ = 1;
    region = 0;
    charge = 0;
    
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.base = [Z, 0.00, 0.00, 0.00, rmsd_3d, occ, region, charge;...
                        Z, 0.50, 0.50, 0.50, rmsd_3d, occ, region, charge];
                    
    atoms = ilc_xtl_build(xtl_parm);

    dz = xtl_parm.c/ncu;
    
    lx = na*xtl_parm.a;
    ly = nb*xtl_parm.b;
    lz = nc*xtl_parm.c;
    
    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end