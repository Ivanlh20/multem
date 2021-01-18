import numpy
import time

import multem
import multem.crystalline_materials


def run():

    ##################### Set system configuration #####################
    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 6
    system_conf.gpu_device = 0

    ################## Load multem default parameter ########$$#########
    input_multem = multem.Input_Multislice(system_conf)  # Load default values

    #################### Set simulation experiment #####################
    # eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52,
    # eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multem.simulation_type = 11

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
    ] = multem.crystalline_materials.Au001_xtl(na, nb, nc, ncu, rmsd_3d)

    ###################### Specimen thickness ##########################
    input_multem.thick_type = (
        2  # eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    )
    input_multem.thick = c / numpy.arange(2, 1000, c)  # Array of thickes (�)

    ###################### x-y sampling ################################
    input_multem.nx = 512
    input_multem.ny = 512
    input_multem.bwl = 0  # Band-width limit, 1: true, 0:false

    #################### Microscope parameters ##########################
    input_multem.E_0 = 300  # Acceleration Voltage (keV)
    input_multem.theta = 0.0  # Till ilumination (�)
    input_multem.phi = 0.0  # Till ilumination (�)

    ###################### Illumination model ##########################
    input_multem.illumination_model = 1  # 1: coherente mode, 4: Numerical integration
    input_multem.temporal_spatial_incoh = (
        1  # 1: Temporal and Spatial, 2: Temporal, 3: Spatial
    )

    ######################## condenser lens ########################
    input_multem.cond_lens_m = 0  # Vortex momentum
    input_multem.cond_lens_c_10 = 14.0312  # Defocus (�)
    input_multem.cond_lens_c_30 = 1e-03  # Third order spherical aberration (mm)
    input_multem.cond_lens_c_50 = 0.00  # Fifth order spherical aberration (mm)
    input_multem.cond_lens_c_12 = 0.0  # Twofold astigmatism (�)
    input_multem.cond_lens_phi_12 = (
        0.0  # Azimuthal angle of the twofold astigmatism (�)
    )
    input_multem.cond_lens_c_23 = 0.0  # Threefold astigmatism (�)
    input_multem.cond_lens_phi_23 = (
        0.0  # Azimuthal angle of the threefold astigmatism (�)
    )
    input_multem.cond_lens_inner_aper_ang = 0.0  # Inner aperture (mrad)
    input_multem.cond_lens_outer_aper_ang = 21.0  # Outer aperture (mrad)

    ######### defocus spread function ############
    dsf_sigma = multem.iehwgd_2_sigma(32)  # from defocus spread to standard deviation
    input_multem.cond_lens_ti_sigma = dsf_sigma  # standard deviation (�)
    input_multem.cond_lens_ti_npts = (
        5  # # of integration points. It will be only used if illumination_model=4
    )

    ########## source spread function ############
    ssf_sigma = multem.hwhm_2_sigma(
        0.45
    )  # half width at half maximum to standard deviation
    input_multem.cond_lens_si_sigma = (
        ssf_sigma  # standard deviation: For parallel ilumination(�^-1) otherwise (�)
    )
    input_multem.cond_lens_si_rad_npts = (
        4  # # of integration points. It will be only used if illumination_model=4
    )

    ######### zero defocus reference ############
    input_multem.cond_lens_zero_defocus_type = (
        "First"
    )  # eZDT_First = 1, eZDT_User_Define = 4
    input_multem.cond_lens_zero_defocus_plane = 0

    ##############################STEM ##############################
    input_multem.scanning_type = 2  # eST_Line = 1, eST_Area = 2
    input_multem.scanning_periodic = (
        1  # 1: true, 0:false (periodic boundary conditions)
    )
    input_multem.scanning_ns = 16
    input_multem.scanning_x0 = 3 * a
    input_multem.scanning_y0 = 3 * b
    input_multem.scanning_xe = 4 * a
    input_multem.scanning_ye = 4 * b

    ############################## Circular Detector ##############################
    input_multem.detector.type = 1  # eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
    input_multem.detector.inner_ang = [40]  # Inner angle(mrad)
    input_multem.detector.outer_ang = [160]  # Outer angle(mrad)

    st = time.perf_counter()
    output_radial_detector = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    ############################## Matrix Detector ##############################
    # This yields the same detector as the radial detector above, to check
    # the consistency of the results

    nxh = input_multem.nx / 2
    nyh = input_multem.ny / 2
    dgx = 1 / input_multem.spec_lx
    dgy = 1 / input_multem.spec_ly

    detector = numpy.zeros(shape=(input_multem.ny, input_multem.nx))
    [gx, gy] = numpy.meshgrid(
        numpy.arange(-nxh, nxh, 1) * dgx, numpy.arange(-nyh, nyh, 1) * dgy
    )
    g = numpy.sqrt(gx ** 2 + gy ** 2)

    g_min = numpy.array(
        [
            multem.mrad_2_rAng(input_multem.E_0, theta)
            for theta in input_multem.detector.inner_ang
        ]
    )
    g_max = numpy.array(
        [
            multem.mrad_2_rAng(input_multem.E_0, theta)
            for theta in input_multem.detector.outer_ang
        ]
    )
    detector[(g_min <= g) & (g < g_max)] = 1

    input_multem.detector.type = 3  # eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
    input_multem.detector.fR = detector

    # Plot Matrix detector

    # imagesc(detector)
    # caxis([0 1])
    # axis image
    # colormap gray

    ########################## run multem ##########################

    st = time.perf_counter()
    output_matrix_detector = multem.tem_simulation(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # ##figure(2) clf
    # for i=1:length(output_radial_detector.data)
    #     suptitle(['Thickness = ' num2str(i)])
    #     sub#plot(1, 3, 1)
    #     imagesc(output_radial_detector.data(i).image_tot(1).image)
    #     title('Radial Detector')
    #     #axis image
    #     #colormap jet
    #     colorbar
    #     sub#plot(1, 3, 2)
    #     imagesc(output_matrix_detector.data(i).image_tot(1).image)
    #     title('Matrix Detector')
    #     #axis image
    #     #colormap jet
    #     colorbar
    #     sub#plot(1, 3, 3)
    #     dI = output_radial_detector.data(i).image_tot(1).image - output_matrix_detector.data(i).image_tot(1).image
    #     dI(dI<1e-6) = 0
    #     imagesc(dI)
    #     title('Difference')
    #     #axis image
    #     #colormap jet
    #     colorbar
    #     pause(0.5)
    # end


if __name__ == "__main__":
    run()
