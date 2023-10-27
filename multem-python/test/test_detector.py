import numpy as np
import pytest
import multem
import pickle


def test_detector():

    detector = multem.Detector()
    grid_1 = multem.Grid_2d()
    grid_1.nx = 10
    grid_1.ny = 20
    grid_2 = multem.Grid_2d()
    grid_2.nx = 11
    grid_2.ny = 21
    detector.type = "Matrix"
    detector.fx = [[1, 2], [3, 4]]
    detector.fR = [[2, 3], [4, 5]]
    detector.grid_1d = [grid_1]
    detector.grid_2d = [grid_2]
    detector.fn = ["Hello", "world"]
    detector.inner_ang = [4, 5, 6, 7]
    detector.outer_ang = [5, 6, 7, 8]

    def check():
        grid_1 = multem.Grid_2d()
        grid_1.nx = 10
        grid_1.ny = 20
        grid_2 = multem.Grid_2d()
        grid_2.nx = 11
        grid_2.ny = 21
        assert detector.type.name == "Matrix"
        assert detector.fx == pytest.approx(np.array([[1, 2], [3, 4]]))
        assert detector.fR == pytest.approx(np.array([[2, 3], [4, 5]]))
        assert len(detector.grid_1d) == 1
        assert len(detector.grid_2d) == 1
        assert detector.grid_1d[0].nx == grid_1.nx
        assert detector.grid_1d[0].ny == grid_1.ny
        assert detector.grid_2d[0].nx == grid_2.nx
        assert detector.grid_2d[0].ny == grid_2.ny
        assert detector.fn == ["Hello", "world"]
        assert detector.inner_ang == pytest.approx(np.array([4, 5, 6, 7]))
        assert detector.outer_ang == pytest.approx(np.array([5, 6, 7, 8]))

    check()

    detector = pickle.loads(pickle.dumps(detector))

    check()
