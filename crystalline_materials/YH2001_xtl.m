function [atoms, lx, ly, lz, a, b, c, dz] = YH2001_xtl(na, nb, nc, ncu, rmsd_3d_Y, rmsd_3d_H)
    xtl_parm.na = na;
    xtl_parm.nb = nb;
    xtl_parm.nc = nc;
    
    xtl_parm.a = 5.2000;
    xtl_parm.b = 5.2000;
    xtl_parm.c = 5.2000;
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;
    
    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];
    
    occ = 1;
    region = 0;
    charge = 0;
    
    %H = 1, Y = 39
    % Z x y z rmsd_3d occupancy charge
    xtl_parm.base = [39, 0.0, 0.0, 0.0, rmsd_3d_Y, occ, region, charge;...
                        39, 0.5, 0.5, 0.0, rmsd_3d_Y, occ, region, charge;...
                        1, 0.25, 0.25, 0.25, rmsd_3d_H, occ, region, charge;...
                        1, 0.75, 0.25, 0.25, rmsd_3d_H, occ, region, charge;...
                        1, 0.25, 0.75, 0.25, rmsd_3d_H, occ, region, charge;...
                        1, 0.75, 0.75, 0.25, rmsd_3d_H, occ, region, charge;...
                        39, 0.0, 0.5, 0.5, rmsd_3d_Y, occ, region, charge;...
                        39, 0.5, 0.0, 0.5, rmsd_3d_Y, occ, region, charge;...
                        1, 0.25, 0.25, 0.75, rmsd_3d_H, occ, region, charge;...
                        1, 0.75, 0.25, 0.75, rmsd_3d_H, occ, region, charge;...
                        1, 0.25, 0.75, 0.75, rmsd_3d_H, occ, region, charge;...
                        1, 0.75, 0.75, 0.75, rmsd_3d_H, occ, region, charge];

    atoms = ilc_xtl_build(xtl_parm);

    dz = xtl_parm.c/ncu;
    lx = na*xtl_parm.a;ly = nb*xtl_parm.b;lz = nc*xtl_parm.c;
    
    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end