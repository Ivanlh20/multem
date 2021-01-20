import pytest
import numpy
import pickle

# from numpy.random import random
import multem

test_data = ["float", "double"]


@pytest.mark.parametrize("dtype", test_data)
def test_output_multislice(dtype):
    output = multem.Output_Multislice(dtype=dtype)
    output.output_type = "image_tot"
    output.ndetector = 1
    output.nx = 2
    output.ny = 3
    output.dx = 1.1
    output.dy = 1.2
    output.dr = 1.3
    output.x = [1, 2, 3]
    output.y = [2, 3, 4]
    output.r = [3, 4, 5]
    output.image_tot = [[[2, 3, 4], [3, 4, 5]], [[4, 5, 6], [8, 8, 9]]]
    output.image_coh = [[[1, 2, 3], [2, 3, 4]], [[3, 4, 5], [6, 7, 8]]]
    output.m2psi_tot = [[8, 9, 10]]
    output.m2psi_coh = [[7, 8, 9]]
    output.psi_coh = [[5 + 5j, 6 + 6j, 7 + 7j]]
    output.V = [[5, 6, 7]]
    output.trans = [[2 + 2j, 3 + 3j, 4 + 4j]]
    output.psi_0 = [[1 + 1j, 2 + 2j, 3 + 3j]]

    def check():
        assert output.output_type.name == "image_tot"
        assert output.ndetector == 1
        assert output.nx == 2
        assert output.ny == 3
        assert output.dx == pytest.approx(1.1)
        assert output.dy == pytest.approx(1.2)
        assert output.dr == pytest.approx(1.3)
        assert output.x == pytest.approx([1, 2, 3])
        assert output.y == pytest.approx([2, 3, 4])
        assert output.r == pytest.approx([3, 4, 5])
        assert output.image_tot[0].image == pytest.approx(
            numpy.array([[2, 3, 4], [3, 4, 5]])
        )
        assert output.image_tot[1].image == pytest.approx(
            numpy.array([[4, 5, 6], [8, 8, 9]])
        )
        assert output.image_coh[0].image == pytest.approx(
            numpy.array([[1, 2, 3], [2, 3, 4]])
        )
        assert output.image_coh[1].image == pytest.approx(
            numpy.array([[3, 4, 5], [6, 7, 8]])
        )
        assert output.m2psi_tot == pytest.approx(numpy.array([[8, 9, 10]]))
        assert output.m2psi_coh == pytest.approx(numpy.array([[7, 8, 9]]))
        assert output.psi_coh == pytest.approx(numpy.array([[5 + 5j, 6 + 6j, 7 + 7j]]))
        assert output.V == pytest.approx(numpy.array([[5, 6, 7]]))
        assert output.trans == pytest.approx(numpy.array([[2 + 2j, 3 + 3j, 4 + 4j]]))
        assert output.psi_0 == pytest.approx(numpy.array([[1 + 1j, 2 + 2j, 3 + 3j]]))

    check()

    output = pickle.loads(pickle.dumps(output))

    check()
