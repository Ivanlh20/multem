import time

import multem
import multem.crystalline_materials


def run():

    lx = 50
    ly = 50
    lz = 20
    Z = 6
    rms_3d = 0.09
    d_min = 1.4
    seed = 1983
    rho = 2.2

    st = time.perf_counter()
    atoms = multem.amorp_spec(lx, ly, lz, d_min, Z, rms_3d, rho, seed)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4),'.r')
    # axis equal

    NA = 6.022140857e23

    d = 0
    lx0 = lx - d
    ly0 = ly - d
    lz0 = lz
    # ii = find((0.5*d<=atoms(:,2))&(atoms(:,2)<=lx-0.5*d)&(0.5*d<=atoms(:,3))&(atoms(:,3)<=ly-0.5*d))
    # atoms = atoms(ii, :)
    density = len(atoms) * 12.011 / (lx0 * ly0 * lz0 * NA * (1e-8) ** 3)

    st = time.perf_counter()
    [r, rdf] = multem.rdf_3d(atoms, 8, 200)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # ##figure(2) clf
    # #plot(r, rdf,'-+r')


if __name__ == "__main__":
    run()
