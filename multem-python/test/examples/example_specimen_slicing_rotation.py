import numpy
import multem
import multem.crystalline_materials


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 4
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.simulation_type = "EWRS"
    input_multem.pn_model = "Still_Atom"
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "Planes"
    input_multem.potential_type = "Lobato_0_12"

    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183
    input_multem.pn_single_conf = 0
    input_multem.pn_nconf = 100

    input_multem.spec_rot_theta = 45
    input_multem.spec_rot_u0 = (1, 0, 0)
    input_multem.spec_rot_center_type = "geometric_center"
    input_multem.spec_rot_center_p = (0, 0, 0)

    na = 8
    nb = 8
    nc = 8
    ncu = 2
    rmsd_3d = 0.085

    [
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        input_multem.spec_lz,
        a,
        b,
        c,
        input_multem.spec_dz,
    ] = multem.crystalline_materials.Au001_xtl(na, nb, nc, ncu, rmsd_3d)

    input_multem.spec_lx = 100
    input_multem.spec_ly = 100
    input_multem.spec_lz = 100

    atoms = numpy.array(input_multem.spec_atoms)
    xc = atoms["x"].max() - atoms["x"].min()
    yc = atoms["y"].max() - atoms["y"].min()
    zc = atoms["z"].max() - atoms["z"].min()
    atoms["x"] += input_multem.spec_lx / 2.0 - xc
    atoms["y"] += input_multem.spec_ly / 2.0 - yc
    atoms["z"] += input_multem.spec_lz / 2.0 - zc
    input_multem.spec_atoms = atoms

    [atoms, Slice] = multem.spec_slicing(input_multem)


if __name__ == "__main__":
    run()
