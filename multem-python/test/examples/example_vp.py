import numpy as np
import time
import multem
import multem.crystalline_materials

from math import log


def run():

    Z = 49
    charge = 0

    Rmin = 1e-02
    Rmax = 5.0
    nR = 512
    dlnR = log(Rmax / Rmin) / (nR - 1)
    R = Rmin * np.exp(np.arange(0, (nR - 1), 1) * dlnR)

    st = time.perf_counter()
    [f1, df1] = multem.vp(1, Z, charge, R)
    [f2, df2] = multem.vp(2, Z, charge, R)
    [f3, df3] = multem.vp(3, Z, charge, R)
    [f4, df4] = multem.vp(4, Z, charge, R)
    [f5, df5] = multem.vp(5, Z, charge, R)
    [f6, df6] = multem.vp(6, Z, charge, R)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
