import time

import multem
import multem.crystalline_materials


def run():

    params = multem.CrystalParameters()
    params.na = 10
    params.nb = 10
    params.nc = 10
    params.a = 4.0780
    params.b = 4.0780
    params.c = 4.0780

    occ = 1
    region = 0
    charge = 0

    rmsd_3d = 0.085

    params.layers = [
        [
            (79, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge),
            (79, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge),
        ],
        [
            (79, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge),
            (79, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge),
        ],
    ]

    st = time.perf_counter()
    atoms = multem.crystal_by_layers(params).spec_atoms
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    na = 1
    nb = 1
    nc = 1
    [atoms, lx, ly, lz, a, b, c, dz] = multem.crystalline_materials.SrTiO3001_xtl(
        na, nb, nc, 2, 0.085
    )


if __name__ == "__main__":
    run()
