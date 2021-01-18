import multem
import multem.crystalline_materials


def run():

    ##################### Set system configuration #####################
    system_conf = multem.System_Configuration()
    system_conf.precision = "float"
    system_conf.device = "device"
    system_conf.cpu_nthread = 1
    system_conf.gpu_device = 0

    input_multem = multem.Input_Multislice(system_conf)  # Load default values

    input_multem.pn_model = (
        3  # ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    )
    input_multem.interaction_model = (
        1  # eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    )
    input_multem.potential_slicing = (
        1  # ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
    )
    input_multem.pn_dim = (1, 1, 1)
    input_multem.pn_seed = 300183
    input_multem.pn_nconf = 3

    input_multem.spec_rot_theta = 0  # final angle
    input_multem.spec_rot_u0 = [1, 0, 0]  # unitary vector
    input_multem.spec_rot_center_type = 1  # 1: geometric center, 2: User define
    input_multem.spec_rot_center_p = [0, 0, 0]  # rotation point

    na = 4
    nb = 4
    nc = 20
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

    # ilm_show_crystal(1, input_multem.spec_atoms)

    input_multem.spec_dz = 5.0

    view
    # get spec slicing
    st = time.perf_counter()
    [atoms, Slice] = multem.spec_slicing(input_multem.toStruct)
    print("Time: %.4f seconds" % (time.perf_counter() - st))
    [natoms, _] = size(atoms)[nslice, _] = size(Slice)

    # for i in range(nslice):

    #    i1 = Slice(i, 5) i2 = Slice(i, 6) ii = numpy.arange(i1, i2, 1)
    #    plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), '.k', atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'or')
    #    #set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1])
    #    title('Atomic positions')
    #    ylabel('y','FontSize',14)
    #    xlabel('x','FontSize',12)
    #    axis equal
    #    i2-i1+1
    #    view([1, 0, 0])
    #    pause(0.1)
    # end

    [size(input_multem.spec_atoms, 1), natoms, nslice]
    [input_multem.spec_lx, input_multem.spec_ly, input_multem.spec_lz]


if __name__ == "__main__":
    run()
