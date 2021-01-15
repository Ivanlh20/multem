import multem
from math import sqrt


def fcc001_xtl(Z, l_c, na, nb, nc, ncu, rmsd_3d):

    params = multem.CrystalParameters()
    params.na = na
    params.nb = nb
    params.nc = nc
    params.a = l_c
    params.b = l_c
    params.c = l_c
    params.alpha = 90
    params.beta = 90
    params.gamma = 90

    occ = 1
    region = 0
    charge = 0

    # Z charge x y z rmsd_3d occupancy region charge
    params.layers = [
        [
            (Z, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge),
            (Z, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge),
        ],
        [
            (Z, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge),
            (Z, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge),
        ],
    ]

    atoms = multem.crystal_by_layers(params)

    dz = params.c / ncu
    lx = na * params.a
    ly = nb * params.b
    lz = nc * params.c

    return atoms, lx, ly, lz, params.a, params.b, params.c, dz


def fcc110_xtl(Z, l_c, na, nb, nc, ncu, rmsd_3d):

    params = multem.CrystalParameters()
    params.na = na
    params.nb = nb
    params.nc = nc
    params.a = l_c / sqrt(2)
    params.b = l_c
    params.c = l_c / sqrt(2)
    params.alpha = 90
    params.beta = 90
    params.gamma = 90

    occ = 1
    region = 0
    charge = 0

    # Z charge x y z rmsd_3d occupancy region charge
    params.layers = [
        [(Z, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge)],
        [(Z, 0.5, 0.5, 0.5, rmsd_3d, occ, region, charge)],
    ]

    atoms = multem.crystal_by_layers(params)

    dz = params.c / ncu
    lx = na * params.a
    ly = nb * params.b
    lz = nc * params.c

    return atoms, lx, ly, lz, params.a, params.b, params.c, dz


def Ag001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/silver/crystal_structure.html

    """
    Z = 47
    a = 4.0853
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Ag110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/silver/crystal_structure.html

    """
    Z = 47
    a = 4.0853
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Al001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/aluminium/crystal_structure.html

    """
    Z = 13
    a = 4.0495
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Al110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/aluminium/crystal_structure.html

    """
    Z = 13
    a = 4.0495
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Au001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/aluminium/crystal_structure.html

    """
    Z = 79
    a = 4.0782
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Au110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/aluminium/crystal_structure.html

    """
    Z = 79
    a = 4.0782
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Cu001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/copper/crystal_structure.html

    """
    Z = 29
    a = 3.6149
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Cu110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/copper/crystal_structure.html

    """
    Z = 29
    a = 3.6149
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Ir001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/iridium/crystal_structure.html

    """
    Z = 77
    a = 3.8390
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Ir110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/iridium/crystal_structure.html

    """
    Z = 77
    a = 3.8390
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Ni001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/nickel/crystal_structure.html

    """
    Z = 28
    a = 3.5240
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Ni110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/nickel/crystal_structure.html

    """
    Z = 28
    a = 3.5240
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Pb001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/lead/crystal_structure.html

    """
    Z = 82
    a = 4.9508
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Pb110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/lead/crystal_structure.html

    """
    Z = 82
    a = 4.9508
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Pd001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/palladium/crystal_structure.html

    """
    Z = 46
    a = 3.8907
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Pd110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/palladium/crystal_structure.html

    """
    Z = 46
    a = 3.8907
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Pt001_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/platinum/crystal_structure.html

    """
    Z = 78
    a = 3.9242
    return fcc001_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)


def Pt110_xtl(na, nb, nc, ncu, rmsd_3d=0.085):
    """
    https://www.webelements.com/platinum/crystal_structure.html

    """
    Z = 78
    a = 3.9242
    return fcc110_xtl(Z, a, na, nb, nc, ncu, rmsd_3d)
