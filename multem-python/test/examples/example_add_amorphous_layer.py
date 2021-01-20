import time

import multem
import multem.crystalline_materials


def run():

    na = 20
    nb = 20
    nc = 30
    ncu = 2
    rmsd_3d = 0.085

    [atoms, lx, ly, _, a, b, c, dz] = multem.crystalline_materials.Au001_xtl(
        na, nb, nc, ncu, rmsd_3d
    )

    lz = 20
    Z = 6
    rms_3d = 0.09
    d_min = 1.4
    seed = 1983
    rho = 2.2
    lay_pos = 1

    st = time.perf_counter()
    atoms = multem.add_amorp_lay(
        atoms, lx, ly, lz, d_min, Z, rms_3d, rho, lay_pos, seed
    )
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
