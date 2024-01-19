"""
Tests for the SetBerendsenPdamp class
"""

import copy
import json
from pathlib import Path

import lmfit
import numpy as np

import pytest
from unittest.mock import MagicMock, patch

from set_berendsen_pdamp import SetBerendsenPdamp

CONFIG = {
    "CORES": 4,
    "TSTART": 1650,
    "SEED": 8607844,
    "PSET": [30000, 0],
    "SIM_TIME_STAGE2": 3,
    "PDAMP_INITIAL": 30000,
    "T_TARGET": 1.0,
    "DT_TOL": 0.001,
    "INDIR": "tests/input",
    "POTENTIAL_FILE": "potential.lmp",
    "LAMMPS_INPUT": {
        "STAGE1": {"TEMPLATE": "stage1_template.lmp", "INPUT": "stage1.lmp"},
        "STAGE2": {"TEMPLATE": "stage2_template.lmp", "INPUT": "stage2.lmp"},
    },
    "OUTDIR": "tests/output",
}

CONFIG_FILE = Path("tests/input/config.json")

with open(CONFIG_FILE, "w", encoding="utf-8") as jf:
    json.dump(CONFIG, jf, indent=4)


@pytest.fixture
def sbp():
    """
    Returns a SetBerendsenPdamp object initialized with CONFIG_FILE.
    """
    return SetBerendsenPdamp(CONFIG_FILE)


def test_init_valid_config(sbp):
    """
    Test that a `SetBerendsenPdamp` object is initialized correctly.
    """

    assert sbp.cores == 4
    assert sbp.pdamp_initial == 30000
    assert sbp.temperature == 1650
    assert sbp.seed == 8607844
    assert sbp.pset == [30000, 0]
    assert sbp.sim_time_stage2 == 3
    assert sbp.t_target == 1.0
    assert sbp.dt_tol == 0.001
    assert sbp.indir == "tests/input"
    assert sbp.potential_file == str(Path("tests/input", "potential.lmp"))
    assert sbp.stage1_template == str(Path("tests/input", "stage1_template.lmp"))
    assert sbp.stage1_input == str(Path("tests/input", "stage1.lmp"))
    assert sbp.stage2_template == str(Path("tests/input", "stage2_template.lmp"))
    assert sbp.stage2_input == str(Path("tests/input", "stage2.lmp"))
    assert sbp.outdir == "tests/output"


def test_init_invalid_config():
    """
    Test that a `SetBerendsenPdamp` object raises a FileNotFoundError if the config file does not exist.
    """

    with pytest.raises(FileNotFoundError):
        SetBerendsenPdamp("invalid_config.json")


def test_init_missing_keys():
    """
    Test that a `SetBerendsenPdamp` object raises a KeyError if the config file is missing a key.
    """

    # Remove a key from the config
    incomplete_config = copy.deepcopy(CONFIG)
    incomplete_config.pop("CORES")
    with open("incomplete_config.json", "w", encoding="utf-8") as jf:
        json.dump(incomplete_config, jf)

    # Test that the SetBerendsenPdamp object raises a KeyError
    with pytest.raises(KeyError):
        SetBerendsenPdamp("incomplete_config.json")
    Path("incomplete_config.json").unlink()


def test_init_no_stage1_keys():
    """
    Test case to verify the behavior of a `SetBerendsenPdamp` object initialization when the
    'STAGE1' key is missing from the config file and the stage1.data file does not exist.
    """

    # Rename stage1.data
    Path("tests/input/stage1.data").rename("tests/input/stage1.data.bak")

    # Remove stage1 keys from the config
    no_stage1_config = copy.deepcopy(CONFIG)
    no_stage1_config["LAMMPS_INPUT"].pop("STAGE1")
    with open("no_stage1_config.json", "w", encoding="utf-8") as jf:
        json.dump(no_stage1_config, jf)

    # Test that the SetBerendsenPdamp object raises an AttributeError
    with pytest.raises(AttributeError):
        SetBerendsenPdamp("no_stage1_config.json")
    Path("no_stage1_config.json").unlink()

    # Rename stage1.data.bak back to stage1.data
    Path("tests/input/stage1.data.bak").rename("tests/input/stage1.data")


# @patch.object(SetBerendsenPdamp, "optimize_pdamp")
# @patch.object(SetBerendsenPdamp, "_plot_fit")
# def test_call(mock_plot_fit, mock_optimize_pdamp, sbp):
#     """
#     Test the __call__ method of the sbp object.
#     """
#     sbp.__call__()
#     mock_optimize_pdamp.assert_called_once()
#     mock_plot_fit.assert_called_once()
def test_call(sbp):
    """
    Test the __call__ method of the sbp object.
    """

    # Create a mock for the SetBerendsenPdamp class
    sbp.optimize_pdamp = MagicMock()
    sbp._save_fit = MagicMock()

    # Call the __call__ method
    sbp()

    # Assert that optimize_pdamp and _save_fit were called
    sbp.optimize_pdamp.assert_called_once()
    sbp._save_fit.assert_called_once()


def test_single_string_replacement(sbp):
    """
    Test case for the `replace_in_template` function which replaces all instances of a single
    string in data with a replacement string.
    """

    data = "Hello [NAME]\nWelcome to [CITY]!"

    replacements = {"[NAME]": "John"}
    outfile = Path("output.txt")
    sbp.replace_in_template(data, replacements, outfile)

    with open(outfile, "r", encoding="utf-8") as f:
        result = "".join(f.readlines())
    outfile.unlink()

    assert result == "Hello John\nWelcome to [CITY]!"


def test_multiple_strings_replacement(sbp):
    """
    Test case for the `replace_in_template` function which replaces all instances of multiple
    strings in data with corresponding replacement strings.
    """

    data = "Hello [NAME]\nWelcome to [CITY]!"
    replacements = {"[NAME]": "John", "[CITY]": "New York"}

    outfile = Path("output.txt")
    sbp.replace_in_template(data, replacements, outfile)

    with open(outfile, "r", encoding="utf-8") as f:
        result = "".join(f.readlines())
    outfile.unlink()

    assert result == "Hello John\nWelcome to New York!"


def test_replace_in_template_empty_data(sbp):
    """
    Test case for the `replace_in_template` function when the data is empty.
    """
    data = ""
    replacements = {"[NAME]": "World"}

    outfile = Path("output.txt")
    sbp.replace_in_template(data, replacements, outfile)

    with open(outfile, "r", encoding="utf-8") as f:
        assert f.read() == ""
    outfile.unlink()


def test_replace_in_template_empty_replacements(sbp):
    """
    Test case to verify the behavior of the `replace_in_template` function when given an empty replacements dictionary.
    """
    data = "Hello, [NAME]!"
    replacements = {}

    outfile = Path("output.txt")
    sbp.replace_in_template(data, replacements, outfile)

    with open(outfile, "r", encoding="utf-8") as f:
        assert f.read() == "Hello, [NAME]!"
    outfile.unlink()


def test_replace_in_template_non_dict_replacements(sbp):
    """
    Test case for the `replace_in_template` function when non-dictionary replacements are provided.
    """
    data = "Hello, [NAME]!"
    replacements = "[NAME]"

    outfile = Path("output.txt")
    with pytest.raises(AttributeError):
        sbp.replace_in_template(data, replacements, outfile)


def test_replace_in_template_non_string_data(sbp):
    """
    Test case for the `replace_in_template` function when non-string data is provided.
    """
    data = 1
    replacements = {"[NAME]": "World"}

    outfile = Path("output.txt")
    with pytest.raises(AttributeError):
        sbp.replace_in_template(data, replacements, outfile)


@patch.object(SetBerendsenPdamp, "replace_in_template")
def test_edit_templates_valid(mock_replace_in_template, sbp):
    """
    Test case for the `edit_templates` function which checks that it calls the `replace_in_template` function twice.
    """
    pdamp = 0.5
    sbp.edit_templates(pdamp)
    assert mock_replace_in_template.call_count == 2


@patch.object(SetBerendsenPdamp, "replace_in_template")
def test_edit_templates_iterable(mock_replace_in_template, sbp):
    """
    Test case for the `edit_templates` function which checks that it calls the
    `replace_in_template` function twice when pdamp is an iterable.
    """
    pdamp = [0.5, 0.6]
    sbp.edit_templates(pdamp)
    assert mock_replace_in_template.call_count == 2


def test_edit_templates_empty_iterable(sbp):
    pdamp = []
    with pytest.raises(ValueError):
        sbp.edit_templates(pdamp)


def test_edit_templates(sbp):
    """
    Test for the `edit_templates` function which compares the edited LAMMPS input files to the expected LAMMPS input files.
    """

    # Call the edit_templates method
    sbp.edit_templates(sbp.pdamp_initial)

    # LAMMPS input file names
    stage1_input = Path(sbp.indir, "stage1.lmp")
    stage2_input = Path(sbp.indir, "stage2.lmp")
    stage1_input_expected = Path(sbp.indir, "stage1_expected.lmp")
    stage2_input_expected = Path(sbp.indir, "stage2_expected.lmp")

    # Read LAMMPS input files
    with open(stage1_input, "r", encoding="utf-8") as f:
        stage1_data = f.read()
    with open(stage2_input, "r", encoding="utf-8") as f:
        stage2_data = f.read()
    with open(stage1_input_expected, "r", encoding="utf-8") as f:
        stage1_data_expected = f.read()
    with open(stage2_input_expected, "r", encoding="utf-8") as f:
        stage2_data_expected = f.read()

    # Check that the LAMMPS input files are correct
    assert stage1_data == stage1_data_expected
    assert stage2_data == stage2_data_expected


@patch("set_berendsen_pdamp.LammpsLibrary")
def test_simulate_valid_stage_number(mock_lammps, sbp):
    """
    Test that the `simulate` function works with a valid stage number (1 or 2).
    """

    sbp.stage1_input = "stage1.in"
    sbp.stage2_input = "stage2.in"
    sbp.data_file = "stage1.data"
    Path(sbp.data_file).touch()  # create the file

    sbp.simulate(1)
    mock_lammps.assert_called_once_with(cores=sbp.cores)
    mock_lammps.return_value.file.assert_called_once_with(sbp.stage1_input)

    mock_lammps.reset_mock()
    sbp.simulate(2)
    mock_lammps.assert_called_once_with(cores=sbp.cores)
    mock_lammps.return_value.file.assert_called_once_with(sbp.stage2_input)


def test_simulate_invalid_stage_number(sbp):
    """
    Test that the `simulate` function raises a ValueError if the stage number is not 1 or 2.
    """

    with pytest.raises(ValueError):
        sbp.simulate(3)


def test_simulate_no_data_file(sbp):
    """
    Test that the `simulate` function raises a FileNotFoundError if stage1.data does not exist after stage 1.
    """

    sbp.stage1_input = "stage1.in"
    sbp.data_file = "stage1.data"
    if Path(sbp.data_file).exists():
        Path(sbp.data_file).unlink()  # delete the file if it exists

    with pytest.raises(FileNotFoundError):
        sbp.simulate(1)


def test_fit_tau_valid(sbp):
    sbp.pressure = np.array([1.0, 0.9, 0.8, 0.7, 0.6])
    sbp.pset = [1.0, 1.0]
    dt = sbp._fit_tau()
    assert isinstance(dt, float)


def test_fit_tau_empty_pressure(sbp):
    sbp.pressure = []
    sbp.pset = [1.0, 1.0]
    with pytest.raises(IndexError):
        sbp._fit_tau()


def test_fit_tau_empty_pset(sbp):
    sbp.pressure = np.array([1.0, 0.9, 0.8, 0.7, 0.6])
    sbp.pset = []
    with pytest.raises(IndexError):
        sbp._fit_tau()


def test_fit_tau_non_iterable_pressure(sbp):
    sbp.pressure = "invalid"
    sbp.pset = [1.0, 1.0]
    with pytest.raises(TypeError):
        sbp._fit_tau()


def test_fit_tau_non_iterable_pset(sbp):
    sbp.pressure = np.array([1.0, 0.9, 0.8, 0.7, 0.6])
    sbp.pset = "invalid"
    with pytest.raises(TypeError):
        sbp._fit_tau()


def test_fit_tau(sbp):
    """
    Test that fitting to real pressure data produce the correct parameters using the `_fit_tau` function.
    Pressure data for testing is in tests/pressure2.dat.
    """

    # Read in pressure data
    pressure_file = Path("tests/pressure2.dat")
    (sbp.time, sbp.pressure) = np.loadtxt(pressure_file, unpack=True)

    # Fit to pressure data
    dt = sbp._fit_tau()

    tau = sbp.fit.params["tau"].value
    t_set = -tau * np.log(0.01)
    p0 = sbp.fit.params["p0"].value

    # Assert that the fit is correct
    assert np.isclose(tau, 0.21718916772531416)
    assert np.isclose(t_set, 1.0001930799281837)
    assert np.isclose(p0, 28333.05245564764)
    assert np.isclose(dt, 0.00019307992818373698)


def test_compute_dt(sbp):
    """
    Test the compute_dt method
    """

    # Mock the methods called in compute_dt
    sbp.edit_templates = MagicMock()
    sbp.simulate = MagicMock()
    sbp._fit_tau = MagicMock(return_value=0.5)

    # Mock the np.loadtxt function to return fixed values
    np.loadtxt = MagicMock(return_value=(np.array([0, 1, 2]), np.array([0, 1, 2])))

    # Run the method with a test value for pdamp
    result = sbp.compute_dt(1.0)

    # Check that the result is as expected
    assert result == 0.5

    # Check that the mocked methods were called with the expected arguments
    sbp.edit_templates.assert_called_once_with(1.0)
    sbp.simulate.assert_called()
    sbp._fit_tau.assert_called_once()

    # Check that np.loadtxt was called with the expected arguments
    np.loadtxt.assert_called_once_with(sbp.pressure_files[1], unpack=True)


def test_residual(sbp):
    """
    Test the _residual method
    """

    # Mock the _pressure_function method to return a constant value
    sbp._pressure_function = lambda params: np.array([1.0, 2.0, 3.0])

    # Set the pressure attribute to a known value
    sbp.pressure = np.array([2.0, 3.0, 4.0])

    # Create a Parameters object
    params = lmfit.Parameters()

    # Call the _residual method
    residuals = sbp._residual(params)

    # Assert that the residuals are as expected
    assert np.array_equal(residuals, np.array([1.0, 1.0, 1.0]))


def test_pressure_function(sbp):
    """
    Test the _pressure_function method
    """

    # Set the time attribute to a known value
    sbp.time = np.array([0.0, 1.0, 2.0])

    # Create a lmfit.Parameters object with known values
    params = lmfit.Parameters()
    params.add("tau", value=1.0)
    params.add("p0", value=2.0)
    params.add("pset", value=3.0)

    # Call the _pressure_function method
    pfit = sbp._pressure_function(params)

    # Assert that the fitted pressure is as expected
    expected_pfit = np.array([2.0, 2.6321205588285577, 2.8646647167633876])
    assert np.allclose(pfit, expected_pfit)


@pytest.mark.filterwarnings("error")
def test_check_f(sbp):
    """
    Test the _check_f method when the f is about 15.
    """

    # Parameters for the fale pressure data
    sbp.pset = [1.0, 1.0]
    tau = 1.0

    # Target relaxation time
    sbp.t_target = -tau * np.log(0.01)

    # Time data
    sbp.time = np.arange(0, 100.1, 0.1)

    def generate_data(f_target):
        """
        Generate some fake pressure data with the desired f value.
        Create a lmfit.Parameters object for "p0".
        """

        p0 = f_target + 1.0
        sbp.pressure = p0 * np.exp(-sbp.time / tau) + sbp.pset[1] * (1.0 - np.exp(-sbp.time / tau))
        p_rand = np.random.normal(0, 1.0, len(sbp.pressure))
        sbp.pressure += p_rand

        params = lmfit.Parameters()
        params.add("p0", value=p0)
        sbp.fit = lambda: None
        sbp.fit.params = params

    # f_target = 15. Assert that no exception was raised
    generate_data(15.0)
    sbp._check_f()
    assert True

    # f_target = 7.0. Assert that a UserWarning was raised
    generate_data(7.0)
    with pytest.raises(UserWarning):
        sbp._check_f()

    # f_target = 2.0. Assert that a ValueError was raised
    generate_data(2.0)
    with pytest.raises(ValueError):
        sbp._check_f()


def test_save_fit(sbp):
    """
    Test the _save_fit method
    """

    # Create a mock for the SetBerendsenPdamp class
    sbp.fit = MagicMock()
    sbp.fit.params = {"tau": MagicMock(value=1.0), "p0": MagicMock(value=2.0)}
    sbp._pressure_function = MagicMock(return_value=[3, 4])
    sbp.temperature = 300
    sbp.pdamp = np.array([[0, 4.0]])
    sbp.pset = [0, 5.0]
    sbp.time = [0, 1]
    sbp.pressure = [6, 7]
    sbp.outdir = "."

    sbp._save_fit()

    outfile = Path(sbp.outdir, "fit.json")

    with open(outfile, "r", encoding="utf-8") as jf:
        json_data = json.load(jf)

    assert json_data == {
        "temperature": 300,
        "pdamp": 4.0,
        "t_set": -1.0 * np.log(0.01),
        "tau": 1.0,
        "P0": 2.0,
        "Pset": 5.0,
        "time": [0, 1],
        "pressure": [6, 7],
        "fit_pressure": [3, 4],
    }

    outfile.unlink()
