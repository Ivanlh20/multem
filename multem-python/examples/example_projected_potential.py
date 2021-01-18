import multem
import multem.crystalline_materials


def run():

    ##################### Set system configuration #####################
    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 4
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)  # Load default values

    input_multem.pn_model = (
        3  # ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    )
    input_multem.interaction_model = (
        1  # eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    )
    input_multem.potential_slicing = (
        2  # ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
    )
    input_multem.potential_type = (
        6
    )  # ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

    input_multem.pn_dim = (1, 1, 1)  # phonon dimensions
    input_multem.pn_seed = 300183  # Random seed(frozen phonon)
    input_multem.pn_single_conf = 1  # 1: true, 0:false
    input_multem.pn_nconf = (
        1  # true: phonon configuration, false: number of frozen phonon configurations
    )

    na = 4
    nb = 4
    nc = 4
    ncu = 2
    rmsd_3d = 0.085

    [
        input_multem.spec_atoms,
        input_multem.spec_lx,
        input_multem.spec_ly,
        input_multem.spec_lz,
        a,
        b,
        c,
        input_multem.spec_dz,
    ] = multem.crystalline_materials.Au001_xtl(na, nb, nc, ncu, rmsd_3d)

    input_multem.nx = 2048
    input_multem.ny = 2048

    # clear multem.spec_slicing
    [atoms, Slice] = multem.spec_slicing(input_multem)

    [natoms, _] = size(atoms)[nslice, _] = size(Slice)
    for i in range(nslice):
        input_multem.islice = islice

        input_multem.system_conf.device = 1  # eD_CPU = 1, eD_GPU = 2
        input_multem.system_conf.precision = "float"  # eP_Float = 1, eP_double = 2
        st = time.perf_counter()
        # clear multem.projected_potential
        ouput_multislice_1 = input_multem.multem.projected_potential
        print("Time: %.4f seconds" % (time.perf_counter() - st))

        input_multem.system_conf.device = "device"  # eD_CPU = 1, eD_GPU = 2
        input_multem.system_conf.precision = "float"  # eP_Float = 1, eP_double = 2
        st = time.perf_counter()
        # clear multem.projected_potential
        ouput_multislice_2 = input_multem.multem.projected_potential
        print("Time: %.4f seconds" % (time.perf_counter() - st))
        mean(abs(ouput_multislice_1.V[:] - ouput_multislice_2.V[:]))

        # ##figure(1)
        # sub#plot(1, 2, 1)
        # imagesc(ouput_multislice_1.V)
        # #colormap gray
        # #axis image
        # sub#plot(1, 2, 2)
        # imagesc(ouput_multislice_2.V)
        # #colormap gray
        # #axis image
        # #disp([min(ouput_multislice_1.V[:]), min(ouput_multislice_2.V[:])])
        # #disp([max(ouput_multislice_1.V[:]), max(ouput_multislice_2.V[:])])
        # pause(0.10)

    end


if __name__ == "__main__":
    run()
