import pytest
import multem
import pickle


def test_grid_2d():

    input = multem.Input()
    scanning = input.scanning
    scanning.type = "Line"
    scanning.pbc = True
    scanning.spxs = False
    scanning.ns = 10
    scanning.x0 = 0.1
    scanning.y0 = 0.2
    scanning.xe = 0.3
    scanning.ye = 0.4

    def check():
        assert input.scanning.type.name == "Line"
        assert input.scanning.pbc == True
        assert input.scanning.spxs == False
        assert input.scanning.ns == 10
        assert input.scanning.x0 == pytest.approx(0.1)
        assert input.scanning.y0 == pytest.approx(0.2)
        assert input.scanning.xe == pytest.approx(0.3)
        assert input.scanning.ye == pytest.approx(0.4)

    check()

    input = pickle.loads(pickle.dumps(input))

    check()

