import multem
import pickle


def test_range_2d():

    range_2d = multem.Range_2d()
    range_2d.ix_0 = 1
    range_2d.ix_e = 2
    range_2d.iy_0 = 3
    range_2d.iy_e = 4
    range_2d.ixy_0 = 5
    range_2d.ixy_e = 6

    def check():
        assert range_2d.ix_0 == 1
        assert range_2d.ix_e == 2
        assert range_2d.iy_0 == 3
        assert range_2d.iy_e == 4
        assert range_2d.ixy_0 == 5
        assert range_2d.ixy_e == 6

    check()

    range_2d = pickle.loads(pickle.dumps(range_2d))

    check()
