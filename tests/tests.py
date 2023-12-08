"""Tests generated with assistance from CodiumAI and GitHub Copilot"""

import json
from pathlib import Path
import sys
from unittest.mock import Mock, patch

import lmfit
import numpy as np
import pytest
from scipy.stats import ks_2samp

from set_berendsen_pdamp import SetBerendsenPdamp

sys.path.append("tests")
from sample_pdamp import multi_sims

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
    return SetBerendsenPdamp(CONFIG_FILE)


def test_init_valid_config(sbp):
    """Test that the SetBerendsenPdamp object is initialized correctly."""

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
    """Test that the SetBerendsenPdamp object raises a FileNotFoundError if the config file does not exist."""

    with pytest.raises(FileNotFoundError):
        obj = SetBerendsenPdamp("invalid_config.json")


def test_init_missing_keys():
    """Test that the SetBerendsenPdamp object raises a KeyError if the config file is missing a key."""

    # Remove a key from the config
    incomplete_config = CONFIG.copy()
    incomplete_config.pop("CORES")
    with open("incomplete_config.json", "w", encoding="utf-8") as jf:
        json.dump(incomplete_config, jf)
    with pytest.raises(KeyError):
        obj = SetBerendsenPdamp("incomplete_config.json")
    Path("incomplete_config.json").unlink()


def test_single_string_replacement(sbp):
    """Replaces all instances of a single string in data with a replacement string and writes to outfile."""

    data = "Hello [NAME]\nWelcome to [CITY]!"
    replacements = {"[NAME]": "John"}
    outfile = Path("output.txt")

    sbp.replace_in_template(data, replacements, outfile)

    with open(outfile, "r", encoding="utf-8") as f:
        result = "".join(f.readlines())
    outfile.unlink()

    assert result == "Hello John\nWelcome to [CITY]!"


def test_multiple_strings_replacement(sbp):
    """Replaces all instances of multiple strings in data with corresponding replacement strings and writes to outfile."""

    data = "Hello [NAME]\nWelcome to [CITY]!"
    replacements = {"[NAME]": "John", "[CITY]": "New York"}
    outfile = Path("output.txt")

    sbp.replace_in_template(data, replacements, outfile)

    with open(outfile, "r", encoding="utf-8") as f:
        result = "".join(f.readlines())
    outfile.unlink()

    assert result == "Hello John\nWelcome to New York!"


def test_edit_templates(sbp):
    """Edit LAMMPS input file templates"""

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
    """Test that the simulate method works with a valid stage number (1 or 2)."""

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
    """The simulate method raises a ValueError if the stage number is not 1 or 2."""

    with pytest.raises(ValueError):
        sbp.simulate(3)


def test_simulate_no_data_file(sbp):
    """The simulate method raises a FileNotFoundError if stage1.data does not exist after stage 1."""

    sbp.stage1_input = "stage1.in"
    sbp.data_file = "stage1.data"
    if Path(sbp.data_file).exists():
        Path(sbp.data_file).unlink()  # delete the file if it exists

    with pytest.raises(FileNotFoundError):
        sbp.simulate(1)


def test_fit_tau(sbp):
    """Test that fitting to pressure data works correctly. Pressure data for testing is in tests/pressure2.dat."""

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


def test_optimization():
    """Since different versions of LAMMPS or a different number of cores might lead to different pdamp values, check that produced pdamp values are from the same distribution as the pre-computed pdamp values in pdamp_samples.json using the Kolmogorov-Smirnov test. Only run stage 2 simulations starting from stage1.data."""

    # Random seeds to use for testing
    seeds = [
        9229241,
        8875157,
        9457236,
        4391786,
        7034636,
        4811723,
        9098824,
        3998610,
        1743382,
        4568358,
    ]

    # Run simulations
    _, data_output = multi_sims(seeds)
    pdamp_test = np.array(data_output["pdamp"])

    # Read in pre-computed pdamp values
    with open(Path("tests/pdamp_samples.json"), "r", encoding="utf-8") as jf:
        pdamp_samples = np.array(json.load(jf)["data"]["pdamp"])

    # Assert that the pdamp values are from the same distribution
    results = ks_2samp(pdamp_test, pdamp_samples)
    assert results.pvalue > 0.05
