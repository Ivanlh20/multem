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

    # eTEMST_EWFS=51, eTEMST_EWRS=52
    input_multem.simulation_type = 52
    input_multem.interaction_model = "Multislice"
    input_multem.potential_slicing = "Planes"
    input_multem.potential_type = "Lobato_0_12"
    input_multem.pn_model = "Frozen_Phonon"

    input_multem.pn_dim = (1, 1, 0)  # phonon dimensions
    input_multem.pn_seed = 300183  # Random seed(frozen phonon)
    input_multem.pn_single_conf = 1  # 1: true, 0:false (extract single configuration)
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
    input_multem.thick = numpy.arange(0, 1000, 2 * c)

    input_multem.nx = 2048
    input_multem.ny = 2048

    ########################### Incident wave ##########################
    input_multem.iw_type = "Plane_Wave"
    input_multem.iw_psi = list(
        test.data.read_psi_0(input_multem.nx, input_multem.ny).flatten()
    )
    input_multem.iw_x = [0.5 * input_multem.spec_lx]  # x position
    input_multem.iw_y = [0.5 * input_multem.spec_ly]  # y position

    ######################## condenser lens ########################
    input_multem.cond_lens_m = 0  # Vortex momentum
    input_multem.cond_lens_c_10 = 1110  # Defocus (�)
    input_multem.cond_lens_c_30 = 3.3  # Third order spherical aberration (mm)
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
    input_multem.cond_lens_outer_aper_ang = 7.50  # Outer aperture (mrad)
    input_multem.cond_lens_ti_sigma = 32  # standard deviation (�)
    input_multem.cond_lens_ti_npts = (
        10  # # of integration points. It will be only used if illumination_model=4
    )
    input_multem.cond_lens_si_sigma = (
        0.2  # standard deviation: For parallel ilumination(�^-1) otherwise (�)
    )
    input_multem.cond_lens_si_rad_npts = (
        8  # # of integration points. It will be only used if illumination_model=4
    )
    input_multem.cond_lens_zero_defocus_type = (
        "First"
    )  # eZDT_First = 1, eZDT_User_Define = 4
    input_multem.cond_lens_zero_defocus_plane = 0

    ######################## Objective lens ########################
    input_multem.obj_lens_m = 0  # Vortex momentum
    input_multem.obj_lens_c_10 = 15.836  # Defocus (�)
    input_multem.obj_lens_c_30 = 1e-03  # Third order spherical aberration (mm)
    input_multem.obj_lens_c_50 = 0.00  # Fifth order spherical aberration (mm)
    input_multem.obj_lens_c_12 = 0.0  # Twofold astigmatism (�)
    input_multem.obj_lens_phi_12 = 0.0  # Azimuthal angle of the twofold astigmatism (�)
    input_multem.obj_lens_c_23 = 0.0  # Threefold astigmatism (�)
    input_multem.obj_lens_phi_23 = (
        0.0  # Azimuthal angle of the threefold astigmatism (�)
    )
    input_multem.obj_lens_inner_aper_ang = 0.0  # Inner aperture (mrad)
    input_multem.obj_lens_outer_aper_ang = 24.0  # Outer aperture (mrad)
    input_multem.obj_lens_ti_sigma = 32  # standard deviation (�)
    input_multem.obj_lens_ti_npts = (
        10  # # of integration points. It will be only used if illumination_model=4
    )
    input_multem.obj_lens_zero_defocus_type = (
        4  # eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
    )
    input_multem.obj_lens_zero_defocus_plane = 5

    input_multem.output_area_ix_0 = 1  # x-starting pixel
    input_multem.output_area_iy_0 = 1  # y-starting pixel
    input_multem.output_area_ix_e = 1  # x-final pixel
    input_multem.output_area_iy_e = 1  # y-final pixel

    # clear multem.wave_function
    st = time.perf_counter()
    ouput_multislice = multem.wave_function(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # ##figure(1)
    # for ithk=1:length(ouput_multislice.thick)
    #     psi_coh = flipud(ouput_multislice.data(ithk).psi_coh)

    #     sub#plot(1, 2, 1)
    #     imagesc(abs(psi_coh)**2)
    #     #colormap gray
    #     #axis image
    #     title(strcat('Intensity, thick = ', num2str(ithk)))
    #     sub#plot(1, 2, 2)
    #     imagesc(angle(psi_coh))
    #     #colormap gray
    #     #axis image
    #     title(strcat('Phase, thick = ', num2str(ithk)))
    #     pause(0.25)
    # end


if __name__ == "__main__":
    run()
