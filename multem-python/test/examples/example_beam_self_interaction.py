import multem
import multem.crystalline_materials


def run():

    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 4
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)

    input_multem.E_0 = 200
    input_multem.theta = 0.00
    input_multem.phi = 0.0

    input_multem.spec_lx = 20
    input_multem.spec_ly = 20

    input_multem.nx = 1024
    input_multem.ny = 1024

    input_multem.iw_type = "Convergent_Wave"
    input_multem.iw_psi = [0]
    input_multem.iw_x = [0.0]
    input_multem.iw_y = [0.0]

    input_multem.cond_lens_m = 0
    input_multem.cond_lens_c_10 = 0
    input_multem.cond_lens_c_30 = 0.002
    input_multem.cond_lens_c_50 = 0.00
    input_multem.cond_lens_c_12 = 0
    input_multem.cond_lens_phi_12 = 0.0
    input_multem.cond_lens_c_23 = 0.0
    input_multem.cond_lens_phi_23 = 0.0
    input_multem.cond_lens_inner_aper_ang = 0.0
    input_multem.cond_lens_outer_aper_ang = 21.0
    input_multem.cond_lens_ti_sigma = 32
    input_multem.cond_lens_ti_npts = 10
    input_multem.cond_lens_si_sigma = 0.2
    input_multem.cond_lens_si_rad_npts = 8

    input_multem.iw_x = [0.5 * input_multem.spec_lx]
    input_multem.iw_y = [0.5 * input_multem.spec_ly]

    df0 = multem.scherzer_defocus(input_multem.E_0, input_multem.cond_lens_c_30)

    input_multem.spec_lx = 20
    input_multem.spec_ly = 20
    input_multem.iw_x = [0.5 * input_multem.spec_lx]
    input_multem.iw_y = [0.5 * input_multem.spec_ly]

    input_multem.cond_lens_c_10 = df0
    output_incident_wave = multem.incident_wave(input_multem)
    psi_i = output_incident_wave.psi_0

    thk = 250
    input_multem.cond_lens_c_10 = df0 + thk
    output_incident_wave = multem.incident_wave(input_multem)
    psi_o = output_incident_wave.psi_0

    input_multem.spec_lx = 50
    input_multem.spec_ly = 50
    input_multem.iw_x = [0.5 * input_multem.spec_lx]
    input_multem.iw_y = [0.5 * input_multem.spec_ly]

    input_multem.cond_lens_c_10 = df0
    output_incident_wave = multem.incident_wave(input_multem)
    psi_i = output_incident_wave.psi_0

    thk = 250
    input_multem.cond_lens_c_10 = df0 + thk
    output_incident_wave = multem.incident_wave(input_multem)
    psi_o = output_incident_wave.psi_0


if __name__ == "__main__":
    run()
