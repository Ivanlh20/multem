import numpy
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

    input_multem.E_0 = 300
    input_multem.theta = 0.0
    input_multem.phi = 0.0

    input_multem.nx = 1024
    input_multem.ny = 1024

    input_multem.spec_lx = 50
    input_multem.spec_ly = 50

    input_multem.iw_type = "Convergent_Wave"
    input_multem.iw_psi = list(
        test.data.read_psi_0(input_multem.nx, input_multem.ny).flatten()
    )
    input_multem.iw_x = [0.5 * input_multem.spec_lx]
    input_multem.iw_y = [0.5 * input_multem.spec_ly]
    input_multem.iw_x = list(10 + numpy.array([0, 25, 0, 25]))
    input_multem.iw_y = list(10 + numpy.array([5, 5, 25, 25]))

    input_multem.cond_lens_m = 0
    input_multem.cond_lens_c_10 = -14.0312
    input_multem.cond_lens_c_30 = 1e-03
    input_multem.cond_lens_c_50 = 0.00
    input_multem.cond_lens_c_12 = 0.0
    input_multem.cond_lens_phi_12 = 0.0
    input_multem.cond_lens_c_23 = 0.0
    input_multem.cond_lens_phi_23 = 0.0
    input_multem.cond_lens_inner_aper_ang = 0
    input_multem.cond_lens_outer_aper_ang = 21.0
    input_multem.cond_lens_ti_sigma = 32
    input_multem.cond_lens_ti_npts = 10
    input_multem.cond_lens_si_sigma = 0.2
    input_multem.cond_lens_si_rad_npts = 8
    input_multem.cond_lens_zero_defocus_type = "First"
    input_multem.cond_lens_zero_defocus_plane = 0
    input_multem.cond_lens_c_10 = multem.scherzer_defocus(
        input_multem.E_0, input_multem.cond_lens_c_30
    )

    input_multem.cond_lens_c_10 = input_multem.cond_lens_c_10

    st = time.perf_counter()
    output_incident_wave = multem.incident_wave(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
