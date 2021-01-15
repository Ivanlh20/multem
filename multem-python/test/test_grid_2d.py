import pytest
import multem
import pickle


def test_grid_2d():

    grid_2d = multem.Grid_2d()
    grid_2d.nx = 10
    grid_2d.ny = 20
    grid_2d.nxh = 30
    grid_2d.nyh = 40
    grid_2d.lx = 100
    grid_2d.ly = 100
    grid_2d.dz = 5
    grid_2d.bwl = True
    grid_2d.pbc_xy = False
    grid_2d.Rx_0 = 0.1
    grid_2d.Ry_0 = 0.2
    grid_2d.dRx = 0.3
    grid_2d.dRy = 0.4
    grid_2d.dgx = 0.5
    grid_2d.dgy = 0.5
    grid_2d.gl2_max = 3.1
    grid_2d.alpha = 20

    def check():
        assert grid_2d.nx == 10
        assert grid_2d.ny == 20
        assert grid_2d.nxh == 30
        assert grid_2d.nyh == 40
        assert grid_2d.lx == pytest.approx(100)
        assert grid_2d.ly == pytest.approx(100)
        assert grid_2d.dz == pytest.approx(5)
        assert grid_2d.bwl == True
        assert grid_2d.pbc_xy == False
        assert grid_2d.Rx_0 == pytest.approx(0.1)
        assert grid_2d.Ry_0 == pytest.approx(0.2)
        assert grid_2d.dRx == pytest.approx(0.3)
        assert grid_2d.dRy == pytest.approx(0.4)
        assert grid_2d.dgx == pytest.approx(0.5)
        assert grid_2d.dgy == pytest.approx(0.5)
        assert grid_2d.gl2_max == pytest.approx(3.1)
        assert grid_2d.alpha == pytest.approx(20)

    check()

    grid_2d = pickle.loads(pickle.dumps(grid_2d))

    check()
