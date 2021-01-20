import numpy
import multem
import multem.crystalline_materials


def run():

    lx = 100
    ly = 100
    lz = 100

    na = 8
    nb = 8
    nc = 8
    ncu = 2
    rmsd_3d = 0.085

    [atoms, _, _, _, _, _, _, _] = multem.crystalline_materials.Au001_xtl(
        na, nb, nc, ncu, rmsd_3d
    )

    atoms = numpy.array(atoms)
    xc = atoms["x"].max() - atoms["x"].min()
    yc = atoms["y"].max() - atoms["y"].min()
    zc = atoms["z"].max() - atoms["z"].min()
    atoms["x"] += lx / 2.0 - xc
    atoms["y"] += ly / 2.0 - yc
    atoms["z"] += lz / 2.0 - zc

    theta = 45
    u0 = (1, 1, 0)
    rot_point_type = "geometric_center"
    p0 = (0, 0, 0)

    atoms_r = multem.spec_rot(atoms, theta, u0, rot_point_type, p0)


if __name__ == "__main__":
    run()
