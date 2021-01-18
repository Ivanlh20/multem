import pytest
import multem.crystalline_materials

test_data = [
    multem.crystalline_materials.Ag001_xtl,
    multem.crystalline_materials.Ag110_xtl,
    multem.crystalline_materials.Al001_xtl,
    multem.crystalline_materials.Al110_xtl,
    multem.crystalline_materials.Au001_xtl,
    multem.crystalline_materials.Au110_xtl,
    multem.crystalline_materials.Cu001_xtl,
    multem.crystalline_materials.Cu110_xtl,
    multem.crystalline_materials.Ir001_xtl,
    multem.crystalline_materials.Ir110_xtl,
    multem.crystalline_materials.Ni001_xtl,
    multem.crystalline_materials.Ni110_xtl,
    multem.crystalline_materials.Pb001_xtl,
    multem.crystalline_materials.Pb110_xtl,
    multem.crystalline_materials.Pd001_xtl,
    multem.crystalline_materials.Pd110_xtl,
    multem.crystalline_materials.Pt001_xtl,
    multem.crystalline_materials.Pt110_xtl,
    multem.crystalline_materials.Si001_xtl,
    multem.crystalline_materials.SrTiO3001_xtl,
]


@pytest.mark.parametrize("Class", test_data)
def test_crystal(Class):
    na = 10
    nb = 10
    nc = 10
    ncu = 1
    atoms, lx, ly, lz, a, b, c, dz = Class(na, nb, nc, ncu)
