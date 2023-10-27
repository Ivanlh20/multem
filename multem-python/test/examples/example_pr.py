import numpy as np
import multem
import multem.crystalline_materials

from math import log


def run():

    Z = 79
    occ = 1
    region = 0
    charge = 0

    rmin = 1e-05
    rmax = 0.1
    nr = 512
    dlnr = log(rmax / rmin) / (nr - 1)
    r = rmin * np.exp(np.arange(0, (nr - 1), 1) * dlnr)

    [f1, df1] = multem.pr(1, Z, charge, r)
    [f2, df2] = multem.pr(2, Z, charge, r)
    [f3, df3] = multem.pr(3, Z, charge, r)
    [f4, df4] = multem.pr(4, Z, charge, r)
    [f5, df5] = multem.pr(5, Z, charge, r)
    [f6, df6] = multem.pr(6, Z, charge, r)


if __name__ == "__main__":
    run()
