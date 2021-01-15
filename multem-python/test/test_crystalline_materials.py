import multem.crystalline_materials


def check_crystal(Func):
    na = 10
    nb = 10
    nc = 10
    ncu = 1
    atoms, lx, ly, lz, a, b, c, dz = Func(na, nb, nc, ncu)


def test_crystalline_materials():
    check_crystal(multem.crystalline_materials.Ag001_xtl)
    check_crystal(multem.crystalline_materials.Ag110_xtl)
    check_crystal(multem.crystalline_materials.Al001_xtl)
    check_crystal(multem.crystalline_materials.Al110_xtl)
    check_crystal(multem.crystalline_materials.Au001_xtl)
    check_crystal(multem.crystalline_materials.Au110_xtl)
    check_crystal(multem.crystalline_materials.Cu001_xtl)
    check_crystal(multem.crystalline_materials.Cu110_xtl)
    check_crystal(multem.crystalline_materials.Ir001_xtl)
    check_crystal(multem.crystalline_materials.Ir110_xtl)
    check_crystal(multem.crystalline_materials.Ni001_xtl)
    check_crystal(multem.crystalline_materials.Ni110_xtl)
    check_crystal(multem.crystalline_materials.Pb001_xtl)
    check_crystal(multem.crystalline_materials.Pb110_xtl)
    check_crystal(multem.crystalline_materials.Pd001_xtl)
    check_crystal(multem.crystalline_materials.Pd110_xtl)
    check_crystal(multem.crystalline_materials.Pt001_xtl)
    check_crystal(multem.crystalline_materials.Pt110_xtl)
