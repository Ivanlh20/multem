function [atoms, lx, ly, lz, a, b, c, dz] = LaAlO3001_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    
    xtl_parm.a = 3.79;
    xtl_parm.b = 3.79;
    xtl_parm.c = 3.79;
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;
    
    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];
    
    occ = 1;
    tag = 0;
    charge = 0;

    % Sr = 38, Ti = 22;O = 8
    % La = 57, Al = 13;O = 8
    % Z x y z rmsd_3d occupancy tag charge
    xtl_parm.base = [57, 0.0, 0.0, 0.0, rmsd_3d, occ, tag, charge;...
                        8, 0.5, 0.5, 0.0, rmsd_3d, occ, tag, charge;...
                        8, 0.0, 0.5, 0.5, rmsd_3d, occ, tag, charge;... 
                        8, 0.5, 0.0, 0.5, rmsd_3d, occ, tag, charge;...
                        13, 0.5, 0.5, 0.5, rmsd_3d, occ, tag, charge];
    
    atoms = ilc_xtl_build(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a;ly = nb*xtl_parm.b;lz = nc*xtl_parm.c;
    
    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end