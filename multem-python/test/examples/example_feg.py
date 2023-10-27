import numpy as np

import multem
import multem.crystalline_materials


def run():

    Z = 50
    occ = 1
    region = 0
    charge = 0

    gmin = 0
    gmax = 12
    ng = 512
    dg = (gmax - gmin) / (ng - 1)
    g = np.arange(gmin, gmax, dg)

    [f1, df1] = multem.feg(1, Z, charge, g)
    [f2, df2] = multem.feg(2, Z, charge, g)
    [f3, df3] = multem.feg(3, Z, charge, g)
    [f4, df4] = multem.feg(4, Z, charge, g)
    [f5, df5] = multem.feg(5, Z, charge, g)
    [f6, df6] = multem.feg(6, Z, charge, g)


if __name__ == "__main__":
    run()
