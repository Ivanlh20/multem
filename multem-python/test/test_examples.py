import importlib.util
import os.path
import pytest

test_data = [
    "examples/example_add_amorphous_layer.py",
    "examples/example_amorphous_specimen.py",
    "examples/example_beam_self_interaction.py",
    "examples/example_check_STEM_input.py",
    "examples/example_crystal_by_layers.py",
    "examples/example_feg.py",
    "examples/example_fxg.py",
    "examples/example_incident_wave_1.py",
    "examples/example_incident_wave_2.py",
    "examples/example_incident_wave_3.py",
    "examples/example_incident_wave_4.py",
    "examples/example_microscope_aberrations.py",
    "examples/example_MULTEM_CBED.py",
    "examples/example_MULTEM_CBEI.py",
    "examples/example_MULTEM_ED.py",
    "examples/example_MULTEM_EELS.py",
    "examples/example_MULTEM_EFTEM.py",
    "examples/example_MULTEM_EWFS.py",
    "examples/example_MULTEM_EWFS_till_illumination.py",
    "examples/example_MULTEM_EWRS.py",
    "examples/example_MULTEM_EWRS_till_illumination.py",
    "examples/example_MULTEM_HCTEM.py",
    "examples/example_MULTEM_HRTEM.py",
    "examples/example_MULTEM_ISTEM.py",
    "examples/example_MULTEM_PED.py",
    "examples/example_MULTEM_STEM_EELS.py",
    "examples/example_MULTEM_STEM_matrix_detector.py",
    "examples/example_MULTEM_STEM.py",
    "examples/example_projected_potential.py",
    "examples/example_propagate.py",
    "examples/example_pr.py",
    "examples/example_specimen_planes.py",
    "examples/example_specimen_rotation.py",
    "examples/example_specimen_slicing_rotation.py",
    "examples/example_specimen_slicing_show.py",
    "examples/example_specimen_slicing_with_amorphous_layer.py",
    "examples/example_specimen_slicing_zone_axis.py",
    "examples/example_specimen_subslicing.py",
    "examples/example_transmission_function.py",
    "examples/example_vp.py",
    "examples/example_vr.py",
    "examples/example_vz.py",
    "examples/example_wave_function.py",
]


def import_module(file_path):
    module_name = "example"
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("file_path", test_data)
def test_example(file_path):

    # Import the module
    package_dir = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(package_dir, file_path)
    module = import_module(file_path)

    # Run the module
    module.run()
