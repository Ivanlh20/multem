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

    input_multem.simulation_type = "PED"

    input_multem.interaction_model = "Multislice"
    input_multem.potential_type = "Lobato_0_12"

    input_multem.potential_slicing = "Planes"

    input_multem.pn_model = "Frozen_Phonon"
    input_multem.pn_coh_contrib = 0
    input_multem.pn_single_conf = 0
    input_multem.pn_nconf = 10
    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183

    na = 16
    nb = 16
    nc = 30
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
    ] = multem.crystalline_materials.Cu001_xtl(na, nb, nc, ncu, rmsd_3d)

    input_multem.thick_type = "Whole_Spec"
    input_multem.thick = np.arange(c, 1000, c)

    input_multem.nx = 1024
    input_multem.ny = 1024
    input_multem.bwl = 0

    input_multem.E_0 = 300
    input_multem.theta = 0.0
    input_multem.phi = 0.0

    input_multem.illumination_model = "Coherent"
    input_multem.temporal_spatial_incoh = "Temporal_Spatial"

    ssf_sigma = multem.mrad_2_sigma(input_multem.E_0, 0.02)

    input_multem.obj_lens_inner_aper_ang = 0.0
    input_multem.obj_lens_outer_aper_ang = 0.0

    input_multem.ped_nrot = 30
    input_multem.ped_theta = 3.0

    st = time.perf_counter()
    output_multislice = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
