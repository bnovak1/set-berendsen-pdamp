"""
Tests generated with assistance from CodiumAI
"""

import json
from pathlib import Path

import numpy as np
import pytest
from pytest_mock import mocker

from set_berendsen_pdamp import SetBerendsenPdamp


class TestSetBerendsenPdamp:
    """
    Test the SetBerendsenPdamp class
    """
    
    def test_instantiation_with_valid_json_config(self):
        """
        Class can be instantiated with a valid JSON config file
        """
    
        # Instantiate the SetBerendsenPdamp class with the valid JSON config file
        config_file = Path("tests/input/config.json")
        sbp = SetBerendsenPdamp(config_file)
    
        # Assert that the class was instantiated successfully
        assert isinstance(sbp, SetBerendsenPdamp)

    def test_single_string_replacement(self):
        """
        Replaces all instances of a single string in data with a replacement string and writes to outfile.
        """

        data = "Hello [NAME]\nWelcome to [CITY]!"
        replacements = {"[NAME]": "John"}
        outfile = Path("output.txt")

        config_file = Path("tests/input/config.json")
        sbp = SetBerendsenPdamp(config_file)
        sbp.replace_in_template(data, replacements, outfile)

        with open(outfile, "r", encoding="utf-8") as f:
            result = ''.join(f.readlines())
        outfile.unlink()
            
        assert result == "Hello John\nWelcome to [CITY]!"

    def test_multiple_strings_replacement(self):
        """
        Replaces all instances of multiple strings in data with corresponding replacement strings and writes to outfile.
        """
        data = "Hello [NAME]\nWelcome to [CITY]!"
        replacements = {"[NAME]": "John", "[CITY]": "New York"}
        outfile = Path("output.txt")

        config_file = Path("tests/input/config.json")
        sbp = SetBerendsenPdamp(config_file)
        sbp.replace_in_template(data, replacements, outfile)

        with open(outfile, "r", encoding="utf-8") as f:
            result = ''.join(f.readlines())
        outfile.unlink()

        assert result == "Hello John\nWelcome to New York!"
      
    def test_edit_templates(self):
        """
        Edit LAMMPS input file templates
        """

        # Instantiate the SetBerendsenPdamp class
        config_file = Path("tests/input/config.json")
        sbp = SetBerendsenPdamp(config_file)
        
        # Value of PDAMP_INITIAL from config.json
        with open(config_file, 'r', encoding="utf-8") as jf:
            pdamp = json.load(jf)["PDAMP_INITIAL"]
        
        # Call the edit_templates method
        sbp.edit_templates(pdamp)
            
        # LAMMPS input file names
        stage1_input = Path("tests/input/stage1.lmp")
        stage2_input = Path("tests/input/stage2.lmp")
        stage1_input_expected = Path("tests/input/stage1_expected.lmp")
        stage2_input_expected = Path("tests/input/stage2_expected.lmp")
        
        # Read LAMMPS input files
        with open(stage1_input, 'r', encoding='utf-8') as f:
            stage1_data = f.read()
        with open(stage2_input, 'r', encoding='utf-8') as f:
            stage2_data = f.read()
        with open(stage1_input_expected, 'r', encoding='utf-8') as f:
            stage1_data_expected = f.read()
        with open(stage2_input_expected, 'r', encoding='utf-8') as f:
            stage2_data_expected = f.read()
            
        # Check that the LAMMPS input files are correct
        assert stage1_data == stage1_data_expected
        assert stage2_data == stage2_data_expected
        
        # Remove the LAMMPS input files that were created
        stage1_input.unlink()
        stage2_input.unlink()
        
    # Test that the class can optimize pdamp and plot the fit without errors
    # def test_optimize_pdamp_and_plot_fit(self, mocker):
    #     # Create a mock for the compute_dt method
    #     mocker.patch.object(SetBerendsenPdamp, 'compute_dt')
    
    #     # Create a mock for the _plot_fit method
    #     mocker.patch.object(SetBerendsenPdamp, '_plot_fit')
    
    #     # Instantiate the SetBerendsenPdamp class
    #     set_berendsen_pdamp = SetBerendsenPdamp("config.json")
    
    #     # Call the __call__ method to optimize pdamp and plot the fit
    #     set_berendsen_pdamp()
    
    #     # Assert that the compute_dt method was called
    #     SetBerendsenPdamp.compute_dt.assert_called()
    
    #     # Assert that the _plot_fit method was called
    #     SetBerendsenPdamp._plot_fit.assert_called()

    # Test that the class can compute dt for a given pdamp without errors
    # def test_compute_dt(self, mocker):
    #     # Create a mock for the edit_templates method
    #     mocker.patch.object(SetBerendsenPdamp, 'edit_templates')
    
    #     # Create a mock for the simulate method
    #     mocker.patch.object(SetBerendsenPdamp, 'simulate')
    
    #     # Create a mock for the _fit_tau method
    #     mocker.patch.object(SetBerendsenPdamp, '_fit_tau')
    
    #     # Create a mock for the np.loadtxt method
    #     mocker.patch.object(np, 'loadtxt')
    
    #     # Instantiate the SetBerendsenPdamp class
    #     set_berendsen_pdamp = SetBerendsenPdamp("config.json")
    
    #     # Call the compute_dt method
    #     dt = set_berendsen_pdamp.compute_dt(0.1)
    
    #     # Assert that the edit_templates method was called
    #     SetBerendsenPdamp.edit_templates.assert_called()
    
    #     # Assert that the simulate method was called
    #     SetBerendsenPdamp.simulate.assert_called()
    
    #     # Assert that the np.loadtxt method was called
    #     np.loadtxt.assert_called()
    
    #     # Assert that the _fit_tau method was called
    #     SetBerendsenPdamp._fit_tau.assert_called()
    
    #     # Assert that the dt value is correct
    #     assert dt == 0.5

