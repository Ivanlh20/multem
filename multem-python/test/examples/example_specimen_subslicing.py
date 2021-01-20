import time
import multem
import multem.crystalline_materials

from math import sqrt, pi


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 1
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.pn_model = "Still_Atom"
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "dz_Sub"
    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183
    input_multem.pn_nconf = 1

    input_multem.spec_rot_theta = 0
    input_multem.spec_rot_u0 = (0, 1, 1)
    input_multem.spec_rot_center_type = "geometric_center"
    input_multem.spec_rot_center_p = (0, 0, 0)

    input_multem.spec_lx = 10
    input_multem.spec_ly = 10
    input_multem.spec_lz = 10
    input_multem.spec_dz = 0.5

    occ = 1
    region = 0
    charge = 0
    input_multem.spec_atoms = [
        (29, 2, 2, 0.0, 0.8, 1.0, 0, charge),
        (29, 6, 2, 0.0, 0.8, 1.0, 0, charge),
    ]
    [
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        lz,
        _,
        _,
        _,
        _,
    ] = multem.crystalline_materials.graphene(1, 1.42, sqrt(0.5 / (8 * pi ** 2)))
    input_multem.spec_dz = 0.5

    st = time.perf_counter()
    input_multem.pn_model = "Still_Atom"
    [atoms0, Slice0] = multem.spec_slicing(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    st = time.perf_counter()
    input_multem.pn_model = "Frozen_Phonon"
    [atoms, Slice] = multem.spec_slicing(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
