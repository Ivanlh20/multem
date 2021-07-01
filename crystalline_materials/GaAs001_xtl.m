function [atoms, lx, ly, lz, a, b, c, dz] = GaAs001_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    
    xtl_parm.a = 5.6537;
    xtl_parm.b = 5.6537;
    xtl_parm.c = 5.6537;
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;
    
    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];
    
    occ = 1;
    tag = 0;
    charge = 0;
    
    % Ga = 31, As = 33;
    % Z x y z rmsd_3d occupancy tag charge
    xtl_parm.base = [31, 0.00, 0.00, 0.00, rmsd_3d, occ, tag, charge;...
                        31, 0.50, 0.50, 0.00, rmsd_3d, occ, tag, charge;...
                        33, 0.25, 0.25, 0.25, rmsd_3d, occ, tag, charge;...
                        33, 0.75, 0.75, 0.25, rmsd_3d, occ, tag, charge;...
                        31, 0.00, 0.50, 0.50, rmsd_3d, occ, tag, charge;...
                        31, 0.50, 0.00, 0.50, rmsd_3d, occ, tag, charge;...
                        33, 0.75, 0.25, 0.75, rmsd_3d, occ, tag, charge;...
                        33, 0.25, 0.75, 0.75, rmsd_3d, occ, tag, charge];

    atoms = ilc_xtl_build(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a;ly = nb*xtl_parm.b;lz = nc*xtl_parm.c;
    
    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end