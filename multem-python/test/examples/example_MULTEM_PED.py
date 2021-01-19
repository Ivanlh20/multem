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
    input_multem = multem.Input_Multislice(system_conf)  # Load default values

    #################### Set simulation experiment #####################
    # eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52,
    # eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multem.simulation_type = 41

    ############## Electron-Specimen interaction model #################
    input_multem.interaction_model = (
        1  # eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    )
    input_multem.potential_type = (
        6
    )  # ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

    ####################### Potential slicing ##########################
    input_multem.potential_slicing = (
        1  # ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
    )

    ############### Electron-Phonon interaction model ##################
    input_multem.pn_model = (
        3  # ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    )
    input_multem.pn_coh_contrib = 0
    input_multem.pn_single_conf = 0  # 1: true, 0:false (extract single configuration)
    input_multem.pn_nconf = (
        10
    )  # true: specific phonon configuration, false: number of frozen phonon configurations
    input_multem.pn_dim = (1, 1, 0)  # phonon dimensions (xyz)
    input_multem.pn_seed = 300183  # Random seed(frozen phonon)

    ####################### Specimen information #######################
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

    ###################### Specimen thickness ##########################
    input_multem.thick_type = (
        1  # eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    )
    input_multem.thick = numpy.arange(c, 1000, c)  # Array of thickes (�)

    ###################### x-y sampling ################################
    input_multem.nx = 1024
    input_multem.ny = 1024
    input_multem.bwl = 0  # Band-width limit, 1: true, 0:false

    #################### Microscope parameters #########################
    input_multem.E_0 = 300  # Acceleration Voltage (keV)
    input_multem.theta = 0.0  # Till ilumination (�)
    input_multem.phi = 0.0  # Till ilumination (�)

    ###################### Illumination model ##########################
    input_multem.illumination_model = (
        1
    )  # 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multem.temporal_spatial_incoh = (
        1  # 1: Temporal and Spatial, 2: Temporal, 3: Spatial
    )

    ######################## condenser lens ########################
    ########## source spread function ############
    ssf_sigma = multem.mrad_2_sigma(
        input_multem.E_0, 0.02
    )  # mrad to standard deviation
    # input_multem.obj_lens_ssf_sigma = (
    #     ssf_sigma  # standard deviation: For parallel ilumination(�^-1) otherwise (�)
    # )
    # input_multem.obj_lens_ssf_npoints = (
    #     4  # # of integration points. It will be only used if illumination_model=4
    # )

    ######################## Objective lens ########################
    ############ aperture radius #################
    input_multem.obj_lens_inner_aper_ang = 0.0  # Inner aperture (mrad)
    input_multem.obj_lens_outer_aper_ang = 0.0  # Outer aperture (mrad)

    ############################### PED ############################
    input_multem.ped_nrot = 30  # number of orientations
    input_multem.ped_theta = 3.0  # Precession angle (degrees)

    st = time.perf_counter()
    output_multislice = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    c = 1e7
    # ##figure(1)
    # for i in range(len(output_multislice.data)):
    #     m2psi_tot = output_multislice.data(i).m2psi_tot
    #     m2psi_tot = log(1+c*m2psi_tot/max(m2psi_tot[:]))

    #     I_min = min(m2psi_tot[:])
    #     I_max = max(m2psi_tot[:])

    #     imagesc(output_multislice.x, output_multislice.y, m2psi_tot, [I_min I_max])
    #     title(strcat('Total intensity -  Thick = ', num2str(i)))
    #     #axis image
    #     #colormap gray
    #     pause(0.5)
    # end


if __name__ == "__main__":
    run()
