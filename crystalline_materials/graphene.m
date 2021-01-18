function [atoms, lx, ly, lz] = graphene(n, a, rms)
    xtl_parm.na = n;
    xtl_parm.nb = round(3*n/sqrt(3));
    xtl_parm.nc = 0;
    xtl_parm.a = 3*a;
    xtl_parm.b = sqrt(3)*a;
    xtl_parm.c = 2;
    xtl_parm.nuLayer = 1;
    occ = 1;
    region = 0;
    charge = 0;
    % C = 6
    % Z x y z rmsd_3d occupancy region charge
    xtl_parm.uLayer(1).atoms = [6, 0.0, 0.0, 0.0, rms, occ, region, charge;...
    6, 1/3, 0.0, 0.0, rms, occ, region, charge; 6, 1/2, 1/2, 0.0, rms, occ, region, charge; 6, 5/6, 1/2, 0.0, rms, occ, region, charge];

    atoms = ilc_crystal_by_lays(xtl_parm);
    lx = xtl_parm.a*xtl_parm.na;
    ly = xtl_parm.b*xtl_parm.nb;
    lz = 0.0;
end