import pytest
import multem
import pickle


def test_grid_2d():

    eels = multem.EELS()
    eels.space = "Real"
    eels.E_0 = 0.1
    eels.E_loss = 1.1
    eels.ge = 2.1
    eels.ge2 = 3.1
    eels.gc = 4.1
    eels.gc2 = 5.1
    eels.m_selection = 2
    eels.collection_angle = 6.1
    eels.channelling_type = "Single_Channelling"
    eels.factor = 7.1
    eels.Z = 2
    eels.x = 8.1
    eels.y = 9.1
    eels.occ = 10.1
    eels.g_collection = 11.1

    def check():
        assert eels.space.name == "Real"
        assert eels.E_0 == pytest.approx(0.1)
        assert eels.E_loss == pytest.approx(1.1)
        assert eels.ge == pytest.approx(2.1)
        assert eels.ge2 == pytest.approx(3.1)
        assert eels.gc == pytest.approx(4.1)
        assert eels.gc2 == pytest.approx(5.1)
        assert eels.m_selection == 2
        assert eels.collection_angle == pytest.approx(6.1)
        assert eels.channelling_type.name == "Single_Channelling"
        assert eels.factor == pytest.approx(7.1)
        assert eels.Z == 2
        assert eels.x == pytest.approx(8.1)
        assert eels.y == pytest.approx(9.1)
        assert eels.occ == pytest.approx(10.1)
        assert eels.g_collection == pytest.approx(11.1)

    check()

    eels = pickle.loads(pickle.dumps(eels))

    check()

