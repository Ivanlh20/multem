import numpy
import time

import multem
import multem.crystalline_materials


def run():

    ##################### Set system configuration #####################
    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 1
    system_conf.gpu_device = 0

    ################## Load multem default parameter ########$$#########
    input_multem = multem.Input_Multislice(system_conf)

    #################### Set simulation experiment #####################
    input_multem.simulation_type = "HRTEM"

    ############## Electron-Specimen interaction model #################
    input_multem.interaction_model = "Multislice"
    input_multem.potential_type = "Lobato_0_12"

    ####################### Potential slicing ##########################
    input_multem.potential_slicing = "Planes"

    ############### Electron-Phonon interaction model ##################
    input_multem.pn_model = "Frozen_Phonon"

    input_multem.pn_coh_contrib = 0
    input_multem.pn_single_conf = 0
    input_multem.pn_nconf = 10
    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183

    ####################### Specimen information #######################
    na = 16
    nb = 16
    nc = 20
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

    ###################### Specimen thickness ##########################
    input_multem.thick_type = "Through_Thick"
    input_multem.thick = numpy.arange(c, 1000, c)

    ###################### x-y sampling ################################
    input_multem.nx = 1024
    input_multem.ny = 1024
    input_multem.bwl = False

    #################### Microscope parameters ##########################
    input_multem.E_0 = 300
    input_multem.theta = 0.0
    input_multem.phi = 0.0

    ###################### Illumination model ##########################
    input_multem.illumination_model = "Coherent"
    input_multem.temporal_spatial_incoh = "Temporal_Spatial"

    ######################## condenser lens ########################
    ########## source spread function ############
    ssf_sigma = multem.mrad_2_sigma(input_multem.E_0, 0.02)
    # input_multem.obj_lens_ssf_sigma = (
    #     ssf_sigma  # standard deviation: For parallel ilumination(�^-1) otherwise (�)
    # )
    # input_multem.obj_lens_ssf_npoints = 4

    ######################## Objective lens ########################
    input_multem.obj_lens_m = 0
    input_multem.obj_lens_c_10 = 20
    input_multem.obj_lens_c_30 = 0.04
    input_multem.obj_lens_c_50 = 0.00
    input_multem.obj_lens_c_12 = 0.0
    input_multem.obj_lens_phi_12 = 0.0
    input_multem.obj_lens_c_23 = 0.0
    input_multem.obj_lens_phi_23 = 0.0
    input_multem.obj_lens_inner_aper_ang = 0.0
    input_multem.obj_lens_outer_aper_ang = 0.0

    ######### defocus spread function ############
    dsf_sigma = multem.iehwgd_2_sigma(32)
    input_multem.obj_lens_ti_sigma = dsf_sigma
    input_multem.obj_lens_ti_npts = 5

    ######### zero defocus reference ############
    input_multem.obj_lens_zero_defocus_type = "First"
    input_multem.obj_lens_zero_defocus_plane = 0

    st = time.perf_counter()
    output_multislice = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # ##figure(1)
    # for i in range(len(output_multislice.data)):
    #     imagesc(output_multislice.data(i).m2psi_tot)
    #     title(strcat('Total intensity -  Thick = ', num2str(output_multislice.thick(i))))
    #     #axis image
    #     #colormap gray
    #     pause(0.25)
    # end


if __name__ == "__main__":
    run()
