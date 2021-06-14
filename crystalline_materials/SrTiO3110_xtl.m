function [atoms, lx, ly, lz, a, b, c, dz] = SrTiO3110_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    
    xtl_parm.a = 3.9050;
    xtl_parm.b = 3.9050*sqrt(2);
    xtl_parm.c = 3.9050*sqrt(2);
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;
    
    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];
    
    occ = 1;
    region = 0;
    charge = 0;
    
    % Sr = 38, Ti = 22;O = 8
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.base = [38, 0.50, 0.00, 0.00, rmsd_3d, occ, region, charge;...
                        8, 0.50, 0.50, 0.00, rmsd_3d, occ, region, charge;...
                        22, 0.0, 0.5, 0.0, rmsd_3d, occ, region, charge;...
                        8, 0.00, 0.25, 0.25, rmsd_3d, occ, region, charge...
                        8, 0.00, 0.75, 0.25, rmsd_3d, occ, region, charge;...
                        8, 0.50, 0.00, 0.50, rmsd_3d, occ, region, charge;...
                        38, 0.50, 0.50, 0.50, rmsd_3d, occ, region, charge;...
                        22, 0.00, 0.00, 0.50, rmsd_3d, occ, region, charge;...
                        8, 0.00, 0.25, 0.75, rmsd_3d, occ, region, charge;...
                        8, 0.00, 0.75, 0.75, rmsd_3d, occ, region, charge];

    atoms = ilc_xtl_build(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a;ly = nb*xtl_parm.b;lz = nc*xtl_parm.c;
    
    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end