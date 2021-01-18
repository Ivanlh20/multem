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
    ####################################################################
    lz = 20
    Z = 6
    rms_3d = 0.09
    d_min = 1.4
    seed = 1983
    rho = 2.2
    lay_pos = 1  # 1: top, 2: bottom

    st = time.perf_counter()
    atoms = multem.add_amorp_lay(
        atoms, lx, ly, lz, d_min, Z, rms_3d, rho, lay_pos, seed
    )
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # ilm_show_crystal(1, atoms)
    # # view([1, 0, 0])
    # zlim([min(atoms[:, 4]), max(atoms[:, 4])])


# disp([lx, ly])
# disp([min(atoms[:, 2]), max(atoms[:, 2])])
# disp([min(atoms[:, 3]), max(atoms[:, 3])])
# disp([min(atoms[:, 4]), max(atoms[:,4])])

# atoms(:, 4) = atoms(:, 4)-min(atoms(:, 4))
# save_atomic_position_pdb('amorphous.pdb', atoms, a, b, c, 90, 90, 90)

if __name__ == "__main__":
    run()
