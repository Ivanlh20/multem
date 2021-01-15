import pytest
import multem
import pickle


def test_grid_2d():

    scanning = multem.Scanning()
    scanning.type = "Line"
    scanning.pbc = True
    scanning.spxs = False
    scanning.ns = 10
    scanning.x0 = 0.1
    scanning.y0 = 0.2
    scanning.xe = 0.3
    scanning.ye = 0.4

    def check():
        assert scanning.type.name == "Line"
        assert scanning.pbc == True
        assert scanning.spxs == False
        assert scanning.ns == 10
        assert scanning.x0 == pytest.approx(0.1)
        assert scanning.y0 == pytest.approx(0.2)
        assert scanning.xe == pytest.approx(0.3)
        assert scanning.ye == pytest.approx(0.4)

    check()

    scanning = pickle.loads(pickle.dumps(scanning))

    check()
