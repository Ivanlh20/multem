import numpy as np
import time

import multem
import multem.crystalline_materials

import test.data


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 4
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.simulation_type = "EWRS"
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "Planes"
    input_multem.potential_type = "Lobato_0_12"
    input_multem.pn_model = "Frozen_Phonon"

    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183
    input_multem.pn_single_conf = True
    input_multem.pn_nconf = 3

    input_multem.bwl = 0
    input_multem.E_0 = 100
    input_multem.theta = 0.01
    input_multem.phi = 0.0

    na = 8
    nb = 8
    nc = 3
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

    input_multem.thick_type = "Whole_Spec"
    input_multem.thick = np.arange(0, 1000, 2 * c)

    input_multem.nx = 2048
    input_multem.ny = 2048

    input_multem.iw_type = "Plane_Wave"
    input_multem.iw_psi = list(
        test.data.read_psi_0(input_multem.nx, input_multem.ny).flatten()
    )
    input_multem.iw_x = [0.5 * input_multem.spec_lx]
    input_multem.iw_y = [0.5 * input_multem.spec_ly]

    input_multem.cond_lens_m = 0
    input_multem.cond_lens_c_10 = 1110
    input_multem.cond_lens_c_30 = 3.3
    input_multem.cond_lens_c_50 = 0.00
    input_multem.cond_lens_c_12 = 0.0
    input_multem.cond_lens_phi_12 = 0.0
    input_multem.cond_lens_c_23 = 0.0
    input_multem.cond_lens_phi_23 = 0.0
    input_multem.cond_lens_inner_aper_ang = 0.0
    input_multem.cond_lens_outer_aper_ang = 7.50
    input_multem.cond_lens_ti_sigma = 32
    input_multem.cond_lens_ti_npts = 10
    input_multem.cond_lens_si_sigma = 0.2
    input_multem.cond_lens_si_rad_npts = 8
    input_multem.cond_lens_zero_defocus_type = "First"
    input_multem.cond_lens_zero_defocus_plane = 0

    input_multem.obj_lens_m = 0
    input_multem.obj_lens_c_10 = 15.836
    input_multem.obj_lens_c_30 = 1e-03
    input_multem.obj_lens_c_50 = 0.00
    input_multem.obj_lens_c_12 = 0.0
    input_multem.obj_lens_phi_12 = 0.0
    input_multem.obj_lens_c_23 = 0.0
    input_multem.obj_lens_phi_23 = 0.0
    input_multem.obj_lens_inner_aper_ang = 0.0
    input_multem.obj_lens_outer_aper_ang = 24.0
    input_multem.obj_lens_ti_sigma = 32
    input_multem.obj_lens_ti_npts = 10
    input_multem.obj_lens_zero_defocus_type = "User_Define"
    input_multem.obj_lens_zero_defocus_plane = 5

    input_multem.output_area_ix_0 = 1
    input_multem.output_area_iy_0 = 1
    input_multem.output_area_ix_e = 1
    input_multem.output_area_iy_e = 1

    st = time.perf_counter()
    ouput_multislice = multem.wave_function(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
