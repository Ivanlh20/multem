import numpy
import time
import multem
import multem.crystalline_materials


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 4
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.pn_model = "Frozen_Phonon"
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "dz_Proj"
    input_multem.potential_type = "Lobato_0_12"

    input_multem.pn_dim = (1, 1, 1)
    input_multem.pn_seed = 300183
    input_multem.pn_single_conf = 1
    input_multem.pn_nconf = 1

    na = 4
    nb = 4
    nc = 4
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

    input_multem.nx = 2048
    input_multem.ny = 2048

    [atoms, Slice] = multem.spec_slicing(input_multem)

    nslice = len(Slice)
    for i in range(nslice):
        input_multem.islice = i

        input_multem.system_conf.device = "host"
        input_multem.system_conf.precision = "float"
        st = time.perf_counter()

        ouput_multislice_1 = multem.projected_potential(input_multem)
        print("Time: %.4f seconds" % (time.perf_counter() - st))

        input_multem.system_conf.device = "device"
        input_multem.system_conf.precision = "float"
        st = time.perf_counter()

        ouput_multislice_2 = multem.projected_potential(input_multem)
        print("Time: %.4f seconds" % (time.perf_counter() - st))
        numpy.mean(
            numpy.abs(
                numpy.array(ouput_multislice_1.V)[:]
                - numpy.array(ouput_multislice_2.V)[:]
            )
        )


if __name__ == "__main__":
    run()
