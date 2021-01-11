import numpy
import pytest
import multem
import pickle


def test_atom_data():

    input = multem.Input()
    atom_data = input.atoms
    atom_data.dz = 5.5
    atom_data.l_x = 100.5
    atom_data.l_y = 200.5
    atom_data.l_z = 300.5
    atom_data.ct_na = 10
    atom_data.ct_nb = 20
    atom_data.ct_nc = 30
    atom_data.ct_a = 10.5
    atom_data.ct_b = 11.5
    atom_data.ct_c = 12.5
    atom_data.ct_x0 = 13.5
    atom_data.ct_y0 = 14.5
    atom_data.amorphous_parameters = [ ( 0.1, 0.2, 0.3), ( 0.4, 0.5, 0.6) ] 
    atom_data.spec_atoms = [ ( 1, 2, 3, 4, 5, 6, 7, 8), ( 2, 3, 4, 5, 6, 7, 8, 9) ]

    def check():
        assert input.atoms.dz == pytest.approx(5.5)
        assert input.atoms.l_x == pytest.approx(100.5)
        assert input.atoms.l_y == pytest.approx(200.5)
        assert input.atoms.l_z == pytest.approx(300.5)
        assert input.atoms.ct_na == 10
        assert input.atoms.ct_nb == 20
        assert input.atoms.ct_nc == 30
        assert input.atoms.ct_a == pytest.approx(10.5)
        assert input.atoms.ct_b == pytest.approx(11.5)
        assert input.atoms.ct_c == pytest.approx(12.5)
        assert input.atoms.ct_x0 == pytest.approx(13.5)
        assert input.atoms.ct_y0 == pytest.approx(14.5)
        assert input.atoms.amorphous_parameters == pytest.approx(
            numpy.array([ ( 0.1, 0.2, 0.3), ( 0.4, 0.5, 0.6) ]))
        assert input.atoms.spec_atoms == pytest.approx(
            numpy.array([ ( 1, 2, 3, 4, 5, 6, 7, 8), ( 2, 3, 4, 5, 6, 7, 8, 9) ]))

    check()

    input = pickle.loads(pickle.dumps(input))

    check()




