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

    input_multem.simulation_type = "EWRS"
    input_multem.pn_model = "Still_Atom"
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "Planes"
    input_multem.potential_type = "Lobato_0_12"

    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183
    input_multem.pn_single_conf = False
    input_multem.pn_nconf = 5

    input_multem.illumination_model = "Partial_Coherent"
    input_multem.temporal_spatial_incoh = "Temporal_Spatial"

    input_multem.bwl = 0

    input_multem.E_0 = 300
    input_multem.theta = 0.0
    input_multem.phi = 0.0

    na = 4
    nb = 4
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
    ] = multem.crystalline_materials.Cu001_xtl(na, nb, nc, ncu, rmsd_3d)

    input_multem.nx = 1024
    input_multem.ny = 1024

    input_multem.iw_type = 4
    input_multem.iw_psi = [0]
    input_multem.iw_x = [input_multem.spec_lx / 2]
    input_multem.iw_y = [input_multem.spec_ly / 2]

    input_multem.cond_lens_si_sigma = 0.2

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
    input_multem.obj_lens_zero_defocus_type = "Last"
    input_multem.obj_lens_zero_defocus_plane = 0

    input_multem.simulation_type = "EWRS"

    st = time.perf_counter()
    output_multislice_0 = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    input_multem.simulation_type = "HRTEM"

    st = time.perf_counter()
    output_multislice_1 = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    input_multem.iw_type = "User_Define_Wave"
    input_multem.iw_psi = output_multislice_0.psi_coh[0]

    st = time.perf_counter()
    output_multislice_2 = multem.microscope_aberrations(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
