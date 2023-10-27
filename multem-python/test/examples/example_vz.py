import numpy as np
import time
import multem
import multem.crystalline_materials

from math import log


def run():

    Z = 79
    charge = 1

    rmin = 1e-02
    rmax = 5.0
    nr = 512
    dlnr = log(rmax / rmin) / (nr - 1)
    r = rmin * np.exp(np.arange(0, (nr - 1), 1) * dlnr)
    z0 = -8.0
    ze = 8.0

    st = time.perf_counter()
    [f1, df1] = multem.vz(1, Z, charge, z0, ze, r)
    [f2, df2] = multem.vz(2, Z, charge, z0, ze, r)
    [f3, df3] = multem.vz(3, Z, charge, z0, ze, r)
    [f4, df4] = multem.vz(4, Z, charge, z0, ze, r)
    [f5, df5] = multem.vz(5, Z, charge, z0, ze, r)
    [f6, df6] = multem.vz(6, Z, charge, z0, ze, r)
    print("Time: %.4f seconds" % (time.perf_counter() - st))


if __name__ == "__main__":
    run()
