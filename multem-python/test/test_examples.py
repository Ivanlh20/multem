import importlib.util
import pytest

test_data = [
    "test.examples.example_add_amorphous_layer",
    "test.examples.example_amorphous_specimen",
    "test.examples.example_beam_self_interaction",
    "test.examples.example_check_STEM_input",
    "test.examples.example_crystal_by_layers",
    "test.examples.example_feg",
    "test.examples.example_fxg",
    "test.examples.example_incident_wave_1",
    "test.examples.example_incident_wave_2",
    "test.examples.example_incident_wave_3",
    "test.examples.example_incident_wave_4",
    "test.examples.example_microscope_aberrations",
    "test.examples.example_MULTEM_CBED",
    "test.examples.example_MULTEM_CBEI",
    "test.examples.example_MULTEM_ED",
    "test.examples.example_MULTEM_EELS",
    "test.examples.example_MULTEM_EFTEM",
    "test.examples.example_MULTEM_EWFS",
    "test.examples.example_MULTEM_EWFS_till_illumination",
    "test.examples.example_MULTEM_EWRS",
    "test.examples.example_MULTEM_EWRS_till_illumination",
    "test.examples.example_MULTEM_HCTEM",
    "test.examples.example_MULTEM_HRTEM",
    "test.examples.example_MULTEM_ISTEM",
    "test.examples.example_MULTEM_PED",
    "test.examples.example_MULTEM_STEM_EELS",
    # "test.examples.example_MULTEM_STEM_matrix_detector",
    "test.examples.example_MULTEM_STEM",
    "test.examples.example_projected_potential",
    "test.examples.example_propagate",
    "test.examples.example_pr",
    "test.examples.example_specimen_planes",
    "test.examples.example_specimen_rotation",
    "test.examples.example_specimen_slicing_rotation",
    "test.examples.example_specimen_slicing_show",
    # "test.examples.example_specimen_slicing_with_amorphous_layer",
    "test.examples.example_specimen_slicing_zone_axis",
    "test.examples.example_specimen_subslicing",
    "test.examples.example_transmission_function",
    "test.examples.example_vp",
    "test.examples.example_vr",
    "test.examples.example_vz",
    "test.examples.example_wave_function",
]


@pytest.mark.parametrize("module_name", test_data)
def test_example(module_name):

    # Import the module
    module = importlib.import_module(module_name)

    # Run the module
    module.run()
