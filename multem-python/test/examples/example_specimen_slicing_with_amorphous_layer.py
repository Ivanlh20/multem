import numpy
import time
import multem
import multem.crystalline_materials


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 12
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.pn_model = "Still_Atom"
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
    ncu = 4
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

    input_multem.spec_dz = a / 2

    lz = 20
    Z = 6
    rms_3d = 0.09
    d_min = 1.4
    seed = 1983
    rho = 2.2
    lay_pos = 2

    z_min = min(numpy.array(input_multem.spec_atoms)["z"])
    z_max = max(numpy.array(input_multem.spec_atoms)["z"])
    st = time.perf_counter()
    input_multem.spec_atoms = multem.add_amorp_lay(
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
        input_multem.spec_amorp(1).z_0 = z_min - lz
        input_multem.spec_amorp(1).z_e = z_min
    else:
        input_multem.spec_amorp(1).z_0 = z_max
        input_multem.spec_amorp(1).z_e = z_max + lz

    input_multem.spec_amorp(1).dz = 2.0

    st = time.perf_counter()
    [atoms, Slice] = multem.spec_slicing(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    st = time.perf_counter()
    [z_planes] = multem.spec_planes(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
