import numpy
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

    input_multem.simulation_type = "STEM"

    input_multem.interaction_model = "Multislice"
    input_multem.potential_type = "Lobato_0_12"

    input_multem.potential_slicing = "Planes"

    input_multem.pn_model = "Frozen_Phonon"
    input_multem.pn_coh_contrib = 0
    input_multem.pn_single_conf = 0
    input_multem.pn_nconf = 10
    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183

    na = 8
    nb = 8
    nc = 5
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

    input_multem.thick_type = "Through_Thick"
    input_multem.thick = c / numpy.arange(2, 1000, c)

    input_multem.nx = 512
    input_multem.ny = 512
    input_multem.bwl = 0

    input_multem.E_0 = 300
    input_multem.theta = 0.0
    input_multem.phi = 0.0

    input_multem.illumination_model = "Coherent"
    input_multem.temporal_spatial_incoh = "Temporal_Spatial"

    input_multem.cond_lens_m = 0
    input_multem.cond_lens_c_10 = 14.0312
    input_multem.cond_lens_c_30 = 1e-03
    input_multem.cond_lens_c_50 = 0.00
    input_multem.cond_lens_c_12 = 0.0
    input_multem.cond_lens_phi_12 = 0.0
    input_multem.cond_lens_c_23 = 0.0
    input_multem.cond_lens_phi_23 = 0.0
    input_multem.cond_lens_inner_aper_ang = 0.0
    input_multem.cond_lens_outer_aper_ang = 21.0

    dsf_sigma = multem.iehwgd_2_sigma(32)
    input_multem.cond_lens_ti_sigma = dsf_sigma
    input_multem.cond_lens_ti_npts = 5

    ssf_sigma = multem.hwhm_2_sigma(0.45)
    input_multem.cond_lens_si_sigma = ssf_sigma
    input_multem.cond_lens_si_rad_npts = 8
    input_multem.cond_lens_si_azm_npts = 12

    input_multem.cond_lens_zero_defocus_type = "First"
    input_multem.cond_lens_zero_defocus_plane = 0

    input_multem.scanning_type = "Area"
    input_multem.scanning_periodic = 1
    input_multem.scanning_square_pxs = 0
    input_multem.scanning_ns = 20
    input_multem.scanning_x0 = 3 * a
    input_multem.scanning_y0 = 3 * b
    input_multem.scanning_xe = 4 * a
    input_multem.scanning_ye = 4 * b

    input_multem.detector.type = 1
    input_multem.detector.inner_ang = [40, 80]
    input_multem.detector.outer_ang = [160, 160]

    st = time.perf_counter()
    output_multislice = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))
    output_multislice


if __name__ == "__main__":
    run()
