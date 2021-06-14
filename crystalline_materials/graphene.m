function [atoms, lx, ly, lz, a, b, c] = graphene(n, da, rms)   
    xtl_parm.a = 3*da;
    xtl_parm.b = sqrt(3)*da;
    xtl_parm.c = 2.0;
    
    xtl_parm.alpha = 90;
    xtl_parm.beta = 90;
    xtl_parm.gamma = 90;

    xtl_parm.na = n;
    xtl_parm.nb = round(n*xtl_parm.a/xtl_parm.b);
    xtl_parm.nc = 0;
    
    xtl_parm.sgn = 1;
    xtl_parm.pbc = false;
    xtl_parm.asym_uc = [];

    occ = 1;
    region = 0;
    charge = 0;
    
    % C = 6
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.base = [6, 0.0, 0.0, 0.0, rms, occ, region, charge;...
                        6, 1/3, 0.0, 0.0, rms, occ, region, charge;...
                        6, 1/2, 1/2, 0.0, rms, occ, region, charge;...
                        6, 5/6, 1/2, 0.0, rms, occ, region, charge];

    atoms = ilc_xtl_build(xtl_parm);
    lx = xtl_parm.a*xtl_parm.na;
    ly = xtl_parm.b*xtl_parm.nb;
    lz = 0.0;

    a = xtl_parm.a;
    b = xtl_parm.b;
    c = xtl_parm.c;
end