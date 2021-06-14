function [atoms, lx, ly, lz, a, b, c, dz] = PdH001_xtl(na, nb, nc, ncu, rmsd_3d)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    
    xtl_parm.a = 4.025;
    xtl_parm.b = 4.025;
    xtl_parm.c = 4.025;
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;
    
    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];
    
    occ = 1;
    region = 0;
    charge = 0;

    % Pd = 46, H = 1
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.base = [46, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge;...
                        46, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge;...
                        1, 0.5, 0.0, 0.0, rmsd_3d, occ, region, charge;...
                        1, 0.0, 0.5, 0.0, rmsd_3d, occ, region, charge;...
                        46, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge;...
                        46, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge;...
                        1, 0.0, 0.0, 0.5, rmsd_3d, occ, region, charge;...
                        1, 0.5, 0.5, 0.5, rmsd_3d, occ, region, charge];
                    
    atoms = ilc_xtl_build(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a;ly = nb*xtl_parm.b;lz = nc*xtl_parm.c;
    
    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end