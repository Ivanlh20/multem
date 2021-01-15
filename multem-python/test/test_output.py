import pytest
import numpy
from numpy.random import random
import multem


def random_image_double(shape=(10, 10)):
    return multem.ImageDouble(
        random(shape[0] * shape[1]).reshape(shape).astype(numpy.float64)
    )


def random_image_complex_double(shape=(10, 10)):
    return multem.ImageComplexDouble(
        random(shape[0] * shape[1]).reshape(shape).astype(numpy.float64)
        + random(shape[0] * shape[1]).reshape(shape).astype(numpy.float64) * 1j
    )


def test_data():

    m2psi_tot = random_image_double()
    m2psi_coh = random_image_double()
    psi_coh = random_image_complex_double()
    image_tot1 = random_image_double()
    image_tot2 = random_image_double()
    image_coh1 = random_image_double()
    image_coh2 = random_image_double()

    d = multem.Data()

    d.m2psi_tot = m2psi_tot
    d.m2psi_coh = m2psi_coh
    d.psi_coh = psi_coh

    d.image_tot = [image_tot1, image_tot2]

    d.image_coh = [image_coh1, image_coh2]

    assert pytest.approx(numpy.array(m2psi_tot)) == numpy.array(d.m2psi_tot)
    assert pytest.approx(numpy.array(m2psi_coh)) == numpy.array(d.m2psi_coh)
    assert pytest.approx(numpy.array(psi_coh)) == numpy.array(d.psi_coh)
    assert pytest.approx(numpy.array(image_tot1)) == numpy.array(d.image_tot[0])
    assert pytest.approx(numpy.array(image_tot2)) == numpy.array(d.image_tot[1])
    assert pytest.approx(numpy.array(image_coh1)) == numpy.array(d.image_coh[0])
    assert pytest.approx(numpy.array(image_coh2)) == numpy.array(d.image_coh[1])


def test_output():

    d1 = multem.Data()
    d2 = multem.Data()
    d3 = multem.Data()
    d4 = multem.Data()
    d5 = multem.Data()
    d1.m2psi_tot = random_image_double()
    d2.m2psi_tot = random_image_double()
    d3.m2psi_tot = random_image_double()
    d4.m2psi_tot = random_image_double()
    d5.m2psi_tot = random_image_double()

    output = multem.Output()
    output.dx = 5
    output.dy = 5
    output.x = [1, 2, 3, 4, 5]
    output.y = [1, 2, 3, 4, 5]
    output.thick = [1, 2, 3, 4, 5]
    output.data = [d1, d2, d3, d4, d5]

    assert pytest.approx(5) == output.dx
    assert pytest.approx(5) == output.dy
    assert pytest.approx([1, 2, 3, 4, 5]) == output.x
    assert pytest.approx([1, 2, 3, 4, 5]) == output.y
    assert pytest.approx([1, 2, 3, 4, 5]) == output.thick
    assert len(output.data) == 5
    assert pytest.approx(numpy.array(d1.m2psi_tot)) == numpy.array(
        output.data[0].m2psi_tot
    )
    assert pytest.approx(numpy.array(d2.m2psi_tot)) == numpy.array(
        output.data[1].m2psi_tot
    )
    assert pytest.approx(numpy.array(d3.m2psi_tot)) == numpy.array(
        output.data[2].m2psi_tot
    )
    assert pytest.approx(numpy.array(d4.m2psi_tot)) == numpy.array(
        output.data[3].m2psi_tot
    )
    assert pytest.approx(numpy.array(d5.m2psi_tot)) == numpy.array(
        output.data[4].m2psi_tot
    )
