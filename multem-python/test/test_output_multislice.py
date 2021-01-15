# import pytest
# import numpy
# from numpy.random import random
import multem


# def random_image_double(shape=(10, 10)):
#     return multem.ImageDouble(
#         random(shape[0] * shape[1]).reshape(shape).astype(numpy.float64)
#     )


# def random_image_complex_double(shape=(10, 10)):
#     return multem.ImageComplexDouble(
#         random(shape[0] * shape[1]).reshape(shape).astype(numpy.float64)
#         + random(shape[0] * shape[1]).reshape(shape).astype(numpy.float64) * 1j
#     )


# def test_data():

#     m2psi_tot = random_image_double()
#     m2psi_coh = random_image_double()
#     psi_coh = random_image_complex_double()
#     image_tot1 = random_image_double()
#     image_tot2 = random_image_double()
#     image_coh1 = random_image_double()
#     image_coh2 = random_image_double()

#     d = multem.Data()

#     d.m2psi_tot = m2psi_tot
#     d.m2psi_coh = m2psi_coh
#     d.psi_coh = psi_coh

#     d.image_tot = [image_tot1, image_tot2]

#     d.image_coh = [image_coh1, image_coh2]

#     assert pytest.approx(numpy.array(m2psi_tot)) == numpy.array(d.m2psi_tot)
#     assert pytest.approx(numpy.array(m2psi_coh)) == numpy.array(d.m2psi_coh)
#     assert pytest.approx(numpy.array(psi_coh)) == numpy.array(d.psi_coh)
#     assert pytest.approx(numpy.array(image_tot1)) == numpy.array(d.image_tot[0])
#     assert pytest.approx(numpy.array(image_tot2)) == numpy.array(d.image_tot[1])
#     assert pytest.approx(numpy.array(image_coh1)) == numpy.array(d.image_coh[0])
#     assert pytest.approx(numpy.array(image_coh2)) == numpy.array(d.image_coh[1])


# def test_output():

#     d1 = multem.Data()
#     d2 = multem.Data()
#     d3 = multem.Data()
#     d4 = multem.Data()
#     d5 = multem.Data()
#     d1.m2psi_tot = random_image_double()
#     d2.m2psi_tot = random_image_double()
#     d3.m2psi_tot = random_image_double()
#     d4.m2psi_tot = random_image_double()
#     d5.m2psi_tot = random_image_double()

#     output = multem.Output()
#     output.dx = 5
#     output.dy = 5
#     output.x = [1, 2, 3, 4, 5]
#     output.y = [1, 2, 3, 4, 5]
#     output.thick = [1, 2, 3, 4, 5]
#     output.data = [d1, d2, d3, d4, d5]

#     assert pytest.approx(5) == output.dx
#     assert pytest.approx(5) == output.dy
#     assert pytest.approx([1, 2, 3, 4, 5]) == output.x
#     assert pytest.approx([1, 2, 3, 4, 5]) == output.y
#     assert pytest.approx([1, 2, 3, 4, 5]) == output.thick
#     assert len(output.data) == 5
#     assert pytest.approx(numpy.array(d1.m2psi_tot)) == numpy.array(
#         output.data[0].m2psi_tot
#     )
#     assert pytest.approx(numpy.array(d2.m2psi_tot)) == numpy.array(
#         output.data[1].m2psi_tot
#     )
#     assert pytest.approx(numpy.array(d3.m2psi_tot)) == numpy.array(
#         output.data[2].m2psi_tot
#     )
#     assert pytest.approx(numpy.array(d4.m2psi_tot)) == numpy.array(
#         output.data[3].m2psi_tot
#     )
#     assert pytest.approx(numpy.array(d5.m2psi_tot)) == numpy.array(
#         output.data[4].m2psi_tot
#     )


def test_output_standard():
    def test(output):

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
        output.image_tot = [[[2, 3, 4], [3, 4, 5]], [4, 5, 6], [8, 8, 9]]
        output.image_coh = [[[1, 2, 3], [2, 3, 4]], [3, 4, 5], [6, 7, 8]]
        output.m2psi_tot = [8, 9, 10]
        output.m2psi_coh = [7, 8, 9]
        output.psi_coh = [5 + 5j, 6 + 6j, 7 + 7j]
        output.V = [5, 6, 7]
        output.trans = [2 + 2j, 3 + 3j, 4 + 4j]
        output.psi_0 = [1 + 1j, 2 + 2j, 3 + 3j]
        output.thk_gpu = [True, False, True]

    # test(multem.OutputF())
    test(multem.OutputD())
