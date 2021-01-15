import pytest
import multem
import pickle


def test_lens():

    lens = multem.Lens()
    lens.m = 1
    lens.c_10 = 1.1
    lens.c_12 = 2.1
    lens.phi_12 = 3.1
    lens.c_21 = 4.1
    lens.phi_21 = 5.1
    lens.c_23 = 6.1
    lens.phi_23 = 7.1
    lens.c_30 = 8.1
    lens.c_32 = 9.1
    lens.phi_32 = 10.1
    lens.c_34 = 11.1
    lens.phi_34 = 12.1
    lens.c_41 = 13.1
    lens.phi_41 = 14.1
    lens.c_43 = 15.1
    lens.phi_43 = 16.1
    lens.c_45 = 17.1
    lens.phi_45 = 18.1
    lens.c_50 = 19.1
    lens.c_52 = 20.1
    lens.phi_52 = 21.1
    lens.c_54 = 22.1
    lens.phi_54 = 23.1
    lens.c_56 = 24.1
    lens.phi_56 = 25.1
    lens.inner_aper_ang = 26.1
    lens.outer_aper_ang = 27.1
    lens.ti_a = 28.1
    lens.ti_sigma = 29.1
    lens.ti_beta = 30.1
    lens.ti_npts = 31
    lens.ti_iehwgd = 32.1
    lens.si_a = 33.1
    lens.si_sigma = 34.1
    lens.si_beta = 35.1
    lens.si_rad_npts = 36
    lens.si_azm_npts = 37
    lens.si_iehwgd = 38.1
    lens.si_theta_c = 39.1
    lens.zero_defocus_type = "User_Define"
    lens.zero_defocus_plane = 40.1
    lens.lambda_ = 41.1

    def check():
        assert lens.m == 1
        assert lens.c_10 == pytest.approx(1.1)
        assert lens.c_12 == pytest.approx(2.1)
        assert lens.phi_12 == pytest.approx(3.1)
        assert lens.c_21 == pytest.approx(4.1)
        assert lens.phi_21 == pytest.approx(5.1)
        assert lens.c_23 == pytest.approx(6.1)
        assert lens.phi_23 == pytest.approx(7.1)
        assert lens.c_30 == pytest.approx(8.1)
        assert lens.c_32 == pytest.approx(9.1)
        assert lens.phi_32 == pytest.approx(10.1)
        assert lens.c_34 == pytest.approx(11.1)
        assert lens.phi_34 == pytest.approx(12.1)
        assert lens.c_41 == pytest.approx(13.1)
        assert lens.phi_41 == pytest.approx(14.1)
        assert lens.c_43 == pytest.approx(15.1)
        assert lens.phi_43 == pytest.approx(16.1)
        assert lens.c_45 == pytest.approx(17.1)
        assert lens.phi_45 == pytest.approx(18.1)
        assert lens.c_50 == pytest.approx(19.1)
        assert lens.c_52 == pytest.approx(20.1)
        assert lens.phi_52 == pytest.approx(21.1)
        assert lens.c_54 == pytest.approx(22.1)
        assert lens.phi_54 == pytest.approx(23.1)
        assert lens.c_56 == pytest.approx(24.1)
        assert lens.phi_56 == pytest.approx(25.1)
        assert lens.inner_aper_ang == pytest.approx(26.1)
        assert lens.outer_aper_ang == pytest.approx(27.1)
        assert lens.ti_a == pytest.approx(28.1)
        assert lens.ti_sigma == pytest.approx(29.1)
        assert lens.ti_beta == pytest.approx(30.1)
        assert lens.ti_npts == 31
        assert lens.ti_iehwgd == pytest.approx(32.1)
        assert lens.si_a == pytest.approx(33.1)
        assert lens.si_sigma == pytest.approx(34.1)
        assert lens.si_beta == pytest.approx(35.1)
        assert lens.si_rad_npts == 36
        assert lens.si_azm_npts == 37
        assert lens.si_iehwgd == pytest.approx(38.1)
        assert lens.si_theta_c == pytest.approx(39.1)
        assert lens.zero_defocus_type.name == "User_Define"
        assert lens.zero_defocus_plane == pytest.approx(40.1)
        assert lens.lambda_ == pytest.approx(41.1)

    check()

    lens = pickle.loads(pickle.dumps(lens))

    check()
