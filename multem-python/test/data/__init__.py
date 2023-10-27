import numpy as np
import os.path
import scipy.io
import scipy.interpolate
from math import pi


def read_m2psi_tot_0():
    return scipy.io.loadmat(os.path.join(os.path.dirname(__file__), "m2psi_tot_0.mat"))[
        "m2psi_tot"
    ]


def read_psi_0(nx, ny):
    psi_0 = read_m2psi_tot_0()
    ny0, nx0 = psi_0.shape
    f = scipy.interpolate.interp2d(np.arange(nx0), np.arange(ny0), psi_0)
    X = np.arange(nx) * (nx0 / nx)
    Y = np.arange(ny) * (ny0 / ny)
    A = f(X, Y)
    P = np.exp(1j * 0.5 * pi * A)
    return A * P
