import pytest
import multem
import numpy
import pickle

def test_input_multislice_standard():
    
    input = multem.Input()

    system_conf = multem.System_Configuration()
    system_conf.device = "device"
    system_conf.precision = "float"
    system_conf.cpu_ncores = 1
    system_conf.cpu_nthread = 2
    system_conf.gpu_device = 3
    system_conf.gpu_nstream = 4

    input.system_conf = system_conf

    input.interaction_model = "Multislice"
    input.potential_type = "Lobato_0_12"
    input.operation_mode = "Normal"
    input.reverse_multislice = False

    input.pn_model = "Still_Atom"
    input.pn_coh_contrib = True
    input.pn_single_conf = False
    input.pn_nconf = 20
    input.pn_dim = (1,1,1)
    input.pn_seed = 40
    
    input.fp_dist = 10
    input.fp_iconf_0 = 11
    input.is_crystal = True

    input.spec_rot_theta = 0
    input.spec_rot_u0 = (1, 2, 3)
    input.spec_rot_center_type = "User_Define"
    input.spec_rot_center_p = (4, 5, 6)

    input.thick_type = "Whole_Spec"
    input.thick = [1.1, 2.2, 3.3, 4.4]

    input.potential_slicing = "Planes"

    input.simulation_type = "HRTEM"

    input.iw_type = "Plane_Wave"
    input.iw_psi = [1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j]
    input.iw_x = [1.1, 2.2, 3.3, 4.4, 5.5]
    input.iw_y = [6.6, 7.7, 8.8, 9.9, 0.0]

    input.E_0 = 3.1
    input.lambda_ = 1
    input.theta = 4.1
    input.phi = 5.9

    input.illumination_model = "Partial_Coherent"
    input.temporal_spatial_incoh = "Temporal_Spatial"

    input.cond_lens.c_10 = 40.5
    input.obj_lens.c_10 = 50.5

    input.eels_fr.E_0 = 10.4

    input.slice_storage = True
    input.mul_sign = 1
    input.Vrl = 0.4
    input.nR = 20
    input.nrot = 30
    input.cdl_var_type = "off"
    input.cdl_var = [ 1, 2, 3 ]
    input.iscan = [ 4, 5, 6 ]
    input.beam_x = [ 7, 8, 9 ]
    input.beam_y = [ 10, 11, 12 ]
    input.islice = 13
    input.dp_Shift = True

    def check():
        assert input.system_conf.device.name == "device"
        assert input.system_conf.precision.name == "float"
        assert input.system_conf.cpu_ncores == 1
        assert input.system_conf.cpu_nthread == 2
        assert input.system_conf.gpu_device == 3
        assert input.system_conf.gpu_nstream == 4

        assert input.interaction_model.name == "Multislice"
        assert input.potential_type.name == "Lobato_0_12"
        assert input.operation_mode.name == "Normal"
        assert input.reverse_multislice == False

        assert input.pn_model.name == "Still_Atom"
        assert input.pn_coh_contrib == True
        assert input.pn_single_conf == False
        assert input.pn_nconf == 20
        assert input.pn_dim == (1,1,1)
        assert input.pn_seed == 40

        assert input.fp_dist == 10
        assert input.fp_iconf_0 == 11
        assert input.is_crystal == True

        assert input.spec_rot_theta == 0
        assert input.spec_rot_u0 == pytest.approx(numpy.array((1, 2, 3)) / numpy.linalg.norm((1,2,3)))
        assert input.spec_rot_center_type.name == "User_Define"
        assert input.spec_rot_center_p == pytest.approx((4, 5, 6))

        assert input.thick_type.name == "Whole_Spec"
        assert input.thick == pytest.approx([1.1, 2.2, 3.3, 4.4])

        assert input.potential_slicing.name == "Planes"

        assert input.simulation_type.name == "HRTEM"

        assert input.iw_type.name == "Plane_Wave"
        assert input.iw_psi == pytest.approx([1 + 1j, 2 + 2j, 3 + 3j, 4 + 4j])
        assert input.iw_x == pytest.approx([1.1, 2.2, 3.3, 4.4, 5.5])
        assert input.iw_y == pytest.approx([6.6, 7.7, 8.8, 9.9, 0.0])

        assert input.E_0 == pytest.approx(3.1)
        assert input.lambda_ == pytest.approx(1)
        assert input.theta == pytest.approx(4.1)
        assert input.phi == pytest.approx(5.9)

        assert input.illumination_model.name == "Partial_Coherent"
        assert input.temporal_spatial_incoh.name == "Temporal_Spatial"

        assert input.cond_lens.c_10 == pytest.approx(40.5)
        assert input.obj_lens.c_10 == pytest.approx(50.5)

        assert input.eels_fr.E_0 == pytest.approx(10.4)

        assert input.slice_storage == True
        assert input.mul_sign == 1
        assert input.Vrl == pytest.approx(0.4)
        assert input.nR == pytest.approx(20)
        assert input.nrot == pytest.approx(30)
        assert input.cdl_var_type.name == "off"
        assert input.cdl_var == pytest.approx([ 1, 2, 3 ])
        assert input.iscan == pytest.approx([ 4, 5, 6 ])
        assert input.beam_x == pytest.approx([ 7, 8, 9 ])
        assert input.beam_y == pytest.approx([ 10, 11, 12 ])
        assert input.islice == 13
        assert input.dp_Shift == True

    check()

    input = pickle.loads(pickle.dumps(input))

    check()


def test_input_multislice_extended():
    pass

    input = multem.Input()

    input.spec_atoms = [(1, 2, 3, 4, 5, 6, 7, 8), (2, 3, 4, 5, 6, 7, 8, 9)]

    input.spec_dz = 50.1
    input.spec_lx = 60.1
    input.spec_ly = 70.1
    input.spec_lz = 80.1
    input.spec_cryst_na = 90
    input.spec_cryst_nb = 100
    input.spec_cryst_nc = 110
    input.spec_cryst_a = 120.1
    input.spec_cryst_b = 130.1
    input.spec_cryst_c = 140.1
    input.spec_cryst_x0 = 150.1
    input.spec_cryst_y0 = 160.1
    input.spec_amorp = [(0.1, 0.2, 0.3), (0.4, 0.5, 0.6)]

    input.nx = 1
    input.ny = 2
    input.bwl = True

    input.cond_lens_m = 1
    input.cond_lens_c_10 = 0.1
    input.cond_lens_c_12 = 0.2
    input.cond_lens_phi_12 = 0.3
    input.cond_lens_c_21 = 0.4
    input.cond_lens_phi_21 = 0.5
    input.cond_lens_c_23 = 0.6
    input.cond_lens_phi_23 = 0.7
    input.cond_lens_c_30 = 0.8
    input.cond_lens_c_32 = 0.9
    input.cond_lens_phi_32 = 1.0
    input.cond_lens_c_34 = 1.1
    input.cond_lens_phi_34 = 1.2
    input.cond_lens_c_41 = 1.3
    input.cond_lens_phi_41 = 1.4
    input.cond_lens_c_43 = 1.5
    input.cond_lens_phi_43 = 1.6
    input.cond_lens_c_45 = 1.7
    input.cond_lens_phi_45 = 1.8
    input.cond_lens_c_50 = 1.9
    input.cond_lens_c_52 = 2.0
    input.cond_lens_phi_52 = 2.1
    input.cond_lens_c_54 = 2.2
    input.cond_lens_phi_54 = 2.3
    input.cond_lens_c_56 = 2.4
    input.cond_lens_phi_56 = 2.5
    input.cond_lens_inner_aper_ang = 2.6
    input.cond_lens_outer_aper_ang = 2.7

    input.cond_lens_ti_a = 0.1 
    input.cond_lens_ti_sigma = 1.1 
    input.cond_lens_ti_beta = 2.1 
    input.cond_lens_ti_npts = 10 
    input.cond_lens_si_a = 4.1 
    input.cond_lens_si_sigma = 6.1 
    input.cond_lens_si_beta = 6.1 
    input.cond_lens_si_rad_npts = 20 
    input.cond_lens_si_azm_npts = 30 

    input.cond_lens_zero_defocus_type = "User_Define"
    input.cond_lens_zero_defocus_plane = 0.123

    input.obj_lens_m = 12
    input.obj_lens_c_10 = 0.1
    input.obj_lens_c_12 = 0.2
    input.obj_lens_phi_12 = 0.3
    input.obj_lens_c_21 = 0.4
    input.obj_lens_phi_21 = 0.5
    input.obj_lens_c_23 = 0.6
    input.obj_lens_phi_23 = 0.7
    input.obj_lens_c_30 = 0.8
    input.obj_lens_c_32 = 0.9
    input.obj_lens_phi_32 = 0.10
    input.obj_lens_c_34 = 0.11
    input.obj_lens_phi_34 = 0.12
    input.obj_lens_c_41 = 0.13
    input.obj_lens_phi_41 = 0.14
    input.obj_lens_c_43 = 0.15
    input.obj_lens_phi_43 = 0.16
    input.obj_lens_c_45 = 0.17
    input.obj_lens_phi_45 = 0.18
    input.obj_lens_c_50 = 0.19
    input.obj_lens_c_52 = 0.20
    input.obj_lens_phi_52 = 0.21
    input.obj_lens_c_54 = 0.22
    input.obj_lens_phi_54 = 0.23
    input.obj_lens_c_56 = 0.24
    input.obj_lens_phi_56 = 0.25
    input.obj_lens_inner_aper_ang = 0.26
    input.obj_lens_outer_aper_ang = 0.27

    input.obj_lens_ti_a = 7.1 
    input.obj_lens_ti_sigma = 6.1 
    input.obj_lens_ti_beta = 8.1 
    input.obj_lens_ti_npts = 13 

    input.obj_lens_zero_defocus_type = "User_Define"
    input.obj_lens_zero_defocus_plane = 1.1

    input.scanning_type = "Line"
    input.scanning_periodic = True
    input.scanning_ns = 20
    input.scanning_x0 = 0.1
    input.scanning_y0 = 0.2
    input.scanning_xe = 0.3
    input.scanning_ye = 0.4

    input.ped_nrot = 10
    input.ped_theta = 4.1

    input.hci_nrot = 10
    input.hci_theta = 4.1

    input.eels_Z = 20
    input.eels_E_loss = 0.9
    input.eels_collection_angle = 10.1
    input.eels_m_selection = 30
    input.eels_channelling_type = "Single_Channelling"

    input.eftem_Z = 20
    input.eftem_E_loss = 0.9
    input.eftem_collection_angle = 10.1
    input.eftem_m_selection = 30
    input.eftem_channelling_type = "Single_Channelling"

    input.output_area_ix_0 = 10
    input.output_area_iy_0 = 20
    input.output_area_ix_e = 30
    input.output_area_iy_e = 40

    def check():
        assert input.spec_atoms == pytest.approx([(1, 2, 3, 4, 5, 6, 7, 8), (2, 3, 4, 5, 6, 7, 8, 9)])

        assert input.spec_dz == pytest.approx(50.1)
        assert input.spec_lx == pytest.approx(60.1)
        assert input.spec_ly == pytest.approx(70.1)
        assert input.spec_lz == pytest.approx(80.1)
        assert input.spec_cryst_na == pytest.approx(90)
        assert input.spec_cryst_nb == pytest.approx(100)
        assert input.spec_cryst_nc == pytest.approx(110)
        assert input.spec_cryst_a == pytest.approx(120.1)
        assert input.spec_cryst_b == pytest.approx(130.1)
        assert input.spec_cryst_c == pytest.approx(140.1)
        assert input.spec_cryst_x0 == pytest.approx(150.1)
        assert input.spec_cryst_y0 == pytest.approx(160.1)
        assert input.spec_amorp == pytest.approx([(0.1, 0.2, 0.3), (0.4, 0.5, 0.6)])

        assert input.nx == 1
        assert input.ny == 2
        assert input.bwl == True

        assert input.cond_lens_m == pytest.approx(1)
        assert input.cond_lens_c_10 == pytest.approx(0.1)
        assert input.cond_lens_c_12 == pytest.approx(0.2)
        assert input.cond_lens_phi_12 == pytest.approx(0.3)
        assert input.cond_lens_c_21 == pytest.approx(0.4)
        assert input.cond_lens_phi_21 == pytest.approx(0.5)
        assert input.cond_lens_c_23 == pytest.approx(0.6)
        assert input.cond_lens_phi_23 == pytest.approx(0.7)
        assert input.cond_lens_c_30 == pytest.approx(0.8)
        assert input.cond_lens_c_32 == pytest.approx(0.9)
        assert input.cond_lens_phi_32 == pytest.approx(1.0)
        assert input.cond_lens_c_34 == pytest.approx(1.1)
        assert input.cond_lens_phi_34 == pytest.approx(1.2)
        assert input.cond_lens_c_41 == pytest.approx(1.3)
        assert input.cond_lens_phi_41 == pytest.approx(1.4)
        assert input.cond_lens_c_43 == pytest.approx(1.5)
        assert input.cond_lens_phi_43 == pytest.approx(1.6)
        assert input.cond_lens_c_45 == pytest.approx(1.7)
        assert input.cond_lens_phi_45 == pytest.approx(1.8)
        assert input.cond_lens_c_50 == pytest.approx(1.9)
        assert input.cond_lens_c_52 == pytest.approx(2.0)
        assert input.cond_lens_phi_52 == pytest.approx(2.1)
        assert input.cond_lens_c_54 == pytest.approx(2.2)
        assert input.cond_lens_phi_54 == pytest.approx(2.3)
        assert input.cond_lens_c_56 == pytest.approx(2.4)
        assert input.cond_lens_phi_56 == pytest.approx(2.5)
        assert input.cond_lens_inner_aper_ang == pytest.approx(2.6)
        assert input.cond_lens_outer_aper_ang == pytest.approx(2.7)

        assert input.cond_lens_ti_a == pytest.approx(0.1 )
        assert input.cond_lens_ti_sigma == pytest.approx(1.1 )
        assert input.cond_lens_ti_beta == pytest.approx(2.1 )
        assert input.cond_lens_ti_npts == pytest.approx(10 )
        assert input.cond_lens_si_a == pytest.approx(4.1 )
        assert input.cond_lens_si_sigma == pytest.approx(6.1 )
        assert input.cond_lens_si_beta == pytest.approx(6.1 )
        assert input.cond_lens_si_rad_npts == pytest.approx(20 )
        assert input.cond_lens_si_azm_npts == pytest.approx(30 )

        assert input.cond_lens_zero_defocus_type.name == pytest.approx("User_Define")
        assert input.cond_lens_zero_defocus_plane == pytest.approx(0.123)

        assert input.obj_lens_m == pytest.approx(12)
        assert input.obj_lens_c_10 == pytest.approx(0.1)
        assert input.obj_lens_c_12 == pytest.approx(0.2)
        assert input.obj_lens_phi_12 == pytest.approx(0.3)
        assert input.obj_lens_c_21 == pytest.approx(0.4)
        assert input.obj_lens_phi_21 == pytest.approx(0.5)
        assert input.obj_lens_c_23 == pytest.approx(0.6)
        assert input.obj_lens_phi_23 == pytest.approx(0.7)
        assert input.obj_lens_c_30 == pytest.approx(0.8)
        assert input.obj_lens_c_32 == pytest.approx(0.9)
        assert input.obj_lens_phi_32 == pytest.approx(0.10)
        assert input.obj_lens_c_34 == pytest.approx(0.11)
        assert input.obj_lens_phi_34 == pytest.approx(0.12)
        assert input.obj_lens_c_41 == pytest.approx(0.13)
        assert input.obj_lens_phi_41 == pytest.approx(0.14)
        assert input.obj_lens_c_43 == pytest.approx(0.15)
        assert input.obj_lens_phi_43 == pytest.approx(0.16)
        assert input.obj_lens_c_45 == pytest.approx(0.17)
        assert input.obj_lens_phi_45 == pytest.approx(0.18)
        assert input.obj_lens_c_50 == pytest.approx(0.19)
        assert input.obj_lens_c_52 == pytest.approx(0.20)
        assert input.obj_lens_phi_52 == pytest.approx(0.21)
        assert input.obj_lens_c_54 == pytest.approx(0.22)
        assert input.obj_lens_phi_54 == pytest.approx(0.23)
        assert input.obj_lens_c_56 == pytest.approx(0.24)
        assert input.obj_lens_phi_56 == pytest.approx(0.25)
        assert input.obj_lens_inner_aper_ang == pytest.approx(0.26)
        assert input.obj_lens_outer_aper_ang == pytest.approx(0.27)

        assert input.obj_lens_ti_a == pytest.approx(7.1 )
        assert input.obj_lens_ti_sigma == pytest.approx(6.1 )
        assert input.obj_lens_ti_beta == pytest.approx(8.1 )
        assert input.obj_lens_ti_npts == pytest.approx(13 )

        assert input.obj_lens_zero_defocus_type.name == "User_Define"
        assert input.obj_lens_zero_defocus_plane == pytest.approx(1.1)

        assert input.scanning_type.name == "Line"
        assert input.scanning_periodic == True
        assert input.scanning_ns == 20
        assert input.scanning_x0 == pytest.approx(0.1)
        assert input.scanning_y0 == pytest.approx(0.2)
        assert input.scanning_xe == pytest.approx(0.3)
        assert input.scanning_ye == pytest.approx(0.4)

        assert input.ped_nrot == 10
        assert input.ped_theta == pytest.approx(4.1)

        assert input.hci_nrot == 10
        assert input.hci_theta == pytest.approx(4.1)

        assert input.eels_Z == 20
        assert input.eels_E_loss == pytest.approx(0.9)
        assert input.eels_collection_angle == pytest.approx(10.1)
        assert input.eels_m_selection == 30
        assert input.eels_channelling_type.name == "Single_Channelling"

        assert input.eftem_Z == 20
        assert input.eftem_E_loss == pytest.approx(0.9)
        assert input.eftem_collection_angle == pytest.approx(10.1)
        assert input.eftem_m_selection == 30
        assert input.eftem_channelling_type.name == "Single_Channelling"

        assert input.output_area_ix_0 == 10
        assert input.output_area_iy_0 == 20
        assert input.output_area_ix_e == 30
        assert input.output_area_iy_e == 40

    check()

    input = pickle.loads(pickle.dumps(input))

    check()

