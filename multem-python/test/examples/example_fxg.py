import numpy

import multem
import multem.crystalline_materials


def run():

    Z = 49
    occ = 1
    region = 0
    charge = 0

    gmin = 0
    gmax = 20
    ng = 512
    dg = (gmax - gmin) / (ng - 1)
    g = numpy.arange(gmin, gmax, dg)

    [f1, df1] = multem.fxg(1, Z, charge, g)
    [f2, df2] = multem.fxg(2, Z, charge, g)
    [f3, df3] = multem.fxg(3, Z, charge, g)
    [f4, df4] = multem.fxg(4, Z, charge, g)
    [f5, df5] = multem.fxg(5, Z, charge, g)
    [f6, df6] = multem.fxg(6, Z, charge, g)


if __name__ == "__main__":
    run()
