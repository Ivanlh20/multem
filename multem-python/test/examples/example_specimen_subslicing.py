import time
import multem
import multem.crystalline_materials

from math import sqrt, pi


def run():

    ##################### Set system configuration #####################
    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 1
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)  # Load default values

    input_multem.pn_model = (
        1
    )  # ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    input_multem.interaction_model = (
        1
    )  # eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    input_multem.potential_slicing = (
        3
    )  # ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
    input_multem.pn_dim = (1, 1, 0)
    input_multem.pn_seed = 300183
    input_multem.pn_nconf = 1

    input_multem.spec_rot_theta = 0  # final angle
    input_multem.spec_rot_u0 = (0, 1, 1)  # unitary vector
    input_multem.spec_rot_center_type = 1  # 1: geometric center, 2: User define
    input_multem.spec_rot_center_p = (0, 0, 0)  # rotation point

    input_multem.spec_lx = 10
    input_multem.spec_ly = 10
    input_multem.spec_lz = 10
    input_multem.spec_dz = 0.5

    occ = 1
    region = 0
    charge = 0
    input_multem.spec_atoms = [
        (29, 2, 2, 0.0, 0.8, 1.0, 0, charge),
        (29, 6, 2, 0.0, 0.8, 1.0, 0, charge),
    ]
    [
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        lz,
        _,
        _,
        _,
        _,
    ] = multem.crystalline_materials.graphene(1, 1.42, sqrt(0.5 / (8 * pi ** 2)))
    input_multem.spec_dz = 0.5

    # get spec slicing
    st = time.perf_counter()
    input_multem.pn_model = 1
    [atoms0, Slice0] = multem.spec_slicing(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # [nslice0, _] = size(Slice0)

    st = time.perf_counter()
    input_multem.pn_model = 3
    [atoms, Slice] = multem.spec_slicing(input_multem)
    print("Time: %.4f seconds" % (time.perf_counter() - st))

    # [nslice, _] = size(Slice)

    # plot(atoms(:, 2), atoms(:, 4), '*k')
    # set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    # title('Atomic positions')
    # ylabel('y','FontSize',14)
    # xlabel('x','FontSize',12)
    # axis equal
    # axis([-2 input_multem.spec_lx+2 -5 input_multem.spec_lz + 5])

    # for i in range(nslice):
    #     hold on
    #     #plot([-2 18], [Slice(i, 1) Slice(i, 1)], '-b', [-2 18], [Slice(i, 2) Slice(i, 2)], '-b')
    #     axis equal
    #     axis([-2 input_multem.spec_lx+2 -5 input_multem.spec_lz + 5])
    # end

    # for i in range(nslice):0
    #     hold on
    #     #plot([-2 input_multem.spec_lx+2], [Slice0(i, 1) Slice0(i, 1)], '-r', [-2 input_multem.spec_lx+2], [Slice0(i, 2) Slice0(i, 2)], '-r')
    #     axis equal
    #     axis([-2 input_multem.spec_lx+2 -5 input_multem.spec_lz + 5])
    # end


if __name__ == "__main__":
    run()
