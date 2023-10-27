import numpy as np
import time
import multem
import multem.crystalline_materials


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 1
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.pn_model = "Frozen_Phonon"
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "Planes"
    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183
    input_multem.pn_nconf = 1

    input_multem.spec_rot_theta = 0
    input_multem.spec_rot_u0 = (1, 0, 0)
    input_multem.spec_rot_center_type = "geometric_center"
    input_multem.spec_rot_center_p = (0, 0, 0)

    na = 6
    nb = 6
    nc = 10
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
    ] = multem.crystalline_materials.Au110_xtl(na, nb, nc, ncu, rmsd_3d)

    [
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        input_multem.spec_lz,
        _,
        _,
        _,
        _,
    ] = multem.crystalline_materials.graphene(10, 2.46, 0.085)

    input_multem.spec_dz = 2

    lz = 20
    Z = 6
    rms_3d = 0.09
    d_min = 1.4
    seed = 1983
    rho = 2.2
    lay_pos = 2

    z_min = np.min(np.array(input_multem.spec_atoms)["z"])
    z_max = np.max(np.array(input_multem.spec_atoms)["z"])
    st = time.perf_counter()
    input_multem.atoms = multem.add_amorp_lay(
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        lz,
        d_min,
        Z,
        rms_3d,
        rho,
        lay_pos,
        seed,
    )
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    if lay_pos == 1:
        input_multem.spec_amorp = [(z_min - lz, z_min, 2.0)]
    else:
        input_multem.spec_amorp = [(z_max, z_max + lz, 2.0)]

    lz = 10
    Z = 6
    rms_3d = 0.09
    d_min = 1.4
    seed = 1983
    rho = 2.2
    lay_pos = 1

    z_min = np.min(np.array(input_multem.spec_atoms)["z"])
    z_max = np.max(np.array(input_multem.spec_atoms)["z"])

    st = time.perf_counter()
    input_multem.atoms = multem.add_amorp_lay(
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        lz,
        d_min,
        Z,
        rms_3d,
        rho,
        lay_pos,
        seed,
    )
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    if lay_pos == 1:
        input_multem.spec_amorp = input_multem.spec_amorp + [(z_min - lz, z_min, 2.0)]
    else:
        input_multem.spec_amorp = input_multem.spec_amorp + [(z_max, z_max + lz, 2.0)]

    st = time.perf_counter()
    z_planes = multem.spec_planes(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
