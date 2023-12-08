"""
Tests generated with assistance from CodiumAI and GitHub Copilot
"""

import json
from pathlib import Path
import sys
from unittest.mock import Mock, patch

import lmfit
import numpy as np
import pytest
from scipy.stats import ks_2samp

from set_berendsen_pdamp import SetBerendsenPdamp
sys.path.append('tests')
from sample_pdamp import multi_sims

CONFIG_FILE = Path("tests/input/config.json")

class TestSetBerendsenPdamp:
    """
    Test the SetBerendsenPdamp class
    """

    def test_instantiation_with_valid_json_config(self):
        """
        Class can be instantiated with a valid JSON config file
        """

        # Instantiate the SetBerendsenPdamp class with the valid JSON config file
        sbp = SetBerendsenPdamp(CONFIG_FILE)

        # Assert that the class was instantiated successfully
        assert isinstance(sbp, SetBerendsenPdamp)

    def test_single_string_replacement(self):
        """
        Replaces all instances of a single string in data with a replacement string and writes to outfile.
        """

        data = "Hello [NAME]\nWelcome to [CITY]!"
        replacements = {"[NAME]": "John"}
        outfile = Path("output.txt")

        sbp = SetBerendsenPdamp(CONFIG_FILE)
        sbp.replace_in_template(data, replacements, outfile)

        with open(outfile, "r", encoding="utf-8") as f:
            result = "".join(f.readlines())
        outfile.unlink()

        assert result == "Hello John\nWelcome to [CITY]!"

    def test_multiple_strings_replacement(self):
        """
        Replaces all instances of multiple strings in data with corresponding replacement strings and writes to outfile.
        """
        data = "Hello [NAME]\nWelcome to [CITY]!"
        replacements = {"[NAME]": "John", "[CITY]": "New York"}
        outfile = Path("output.txt")

        sbp = SetBerendsenPdamp(CONFIG_FILE)
        sbp.replace_in_template(data, replacements, outfile)

        with open(outfile, "r", encoding="utf-8") as f:
            result = "".join(f.readlines())
        outfile.unlink()

        assert result == "Hello John\nWelcome to New York!"

    def test_edit_templates(self):
        """
        Edit LAMMPS input file templates
        """
        
        # Make sure that the random seed is correct in the config file
        with open(CONFIG_FILE, "r", encoding="utf-8") as jf:
            config = json.load(jf)
        config["SEED"] = 8607844
        with open(CONFIG_FILE, "w", encoding="utf-8") as jf:
            json.dump(config, jf, indent=4)

        # Instantiate the SetBerendsenPdamp class
        sbp = SetBerendsenPdamp(CONFIG_FILE)

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

@patch('set_berendsen_pdamp.LammpsLibrary')
def test_simulate_valid_stage_number(mock_lammps):
    """
    Test that the simulate method works with a valid stage number (1 or 2).
    """
    
    sbp = SetBerendsenPdamp(CONFIG_FILE)
    sbp.stage1_input = 'stage1.in'
    sbp.stage2_input = 'stage2.in'
    sbp.data_file = 'stage1.data'
    Path(sbp.data_file).touch()  # create the file

    sbp.simulate(1)
    mock_lammps.assert_called_once_with(cores=sbp.cores)
    mock_lammps.return_value.file.assert_called_once_with(sbp.stage1_input)

    mock_lammps.reset_mock()
    sbp.simulate(2)
    mock_lammps.assert_called_once_with(cores=sbp.cores)
    mock_lammps.return_value.file.assert_called_once_with(sbp.stage2_input)

def test_simulate_invalid_stage_number():
    """
    The simulate method raises a ValueError if the stage number is not 1 or 2.
    """
    sbp = SetBerendsenPdamp(CONFIG_FILE)
    with pytest.raises(ValueError):
        sbp.simulate(3)

@patch('set_berendsen_pdamp.LammpsLibrary')
def test_simulate_no_data_file(mock_lammps):
    """
    The simulate method raises a FileNotFoundError if stage1.data does not exist after stage 1.
    """
    sbp = SetBerendsenPdamp(CONFIG_FILE)
    sbp.stage1_input = 'stage1.in'
    sbp.data_file = 'stage1.data'
    if Path(sbp.data_file).exists():
        Path(sbp.data_file).unlink()  # delete the file if it exists

    with pytest.raises(FileNotFoundError):
        sbp.simulate(1)
        
    def test_fit_tau(self):
        """
        Test that fitting to pressure data works correctly. Pressure data for testing is in tests/pressure2.dat.
        """

        # Instantiate the SetBerendsenPdamp class
        sbp = SetBerendsenPdamp(CONFIG_FILE)
        
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

    def test_optimization(self):
        """
        Since different versions of LAMMPS or a different number of cores might lead to different pdamp values, check that produced pdamp values are from the same distribution as the pre-computed pdamp values in pdamp_samples.json using the Kolmogorov-Smirnov test. Only run stage 2 simulations starting from stage1.data.
        """

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