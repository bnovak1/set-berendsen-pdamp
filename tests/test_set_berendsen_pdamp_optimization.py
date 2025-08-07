"""
Test the actual optimization of pdamp.
Since different versions of LAMMPS or a different number of cores might lead to different pdamp values,
check that produced pdamp values are from the same distribution as the pre-computed pdamp
values in pdamp_samples.json using the Kolmogorov-Smirnov test.
Only run stage 2 simulations starting from stage1.data.
"""

import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.stats import ks_2samp

sys.path.append("tests")
from sample_pdamp import multi_sims

@pytest.fixture
def optimization_config():
    """Create a temporary optimization config file for testing."""
    config = {
        "CORES": 2,
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

    # Write config file
    config_file = Path("tests/input/optimization_config.json")
    with open(config_file, "w", encoding="utf-8") as jf:
        json.dump(config, jf, indent=4)

    # Remove the leading "tests/" from the path since it will get added back in multi_sims
    config_file = config_file.relative_to(Path("tests"))
    
    yield config_file
    
    # Cleanup
    config_file = Path("tests/input/optimization_config.json")
    if config_file.exists():
        config_file.unlink()

@pytest.fixture
def seeds():
    """Random seeds to use for testing."""
    random_seeds = [
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
    
    yield random_seeds
    
    # Cleanup
    random_seeds.clear()
    del random_seeds

def test_optimization(optimization_config, seeds):
    """
    Test the actual optimization of pdamp.
    Since different versions of LAMMPS or a different number of cores might lead to different pdamp values,
    check that produced pdamp values are from the same distribution as the pre-computed pdamp
    values in pdamp_samples.json using the Kolmogorov-Smirnov test.
    Only run stage 2 simulations starting from stage1.data.
    """

    # Run simulations
    _, data_output = multi_sims(seeds, infile=optimization_config)
    pdamp_test = np.array(data_output["pdamp"])

    # Read in pre-computed pdamp values
    with open(Path("tests/pdamp_samples.json"), "r", encoding="utf-8") as jf:
        pdamp_samples = np.array(json.load(jf)["data"]["pdamp"])

    # Assert that the pdamp values are from the same distribution
    results = ks_2samp(pdamp_test, pdamp_samples)
    assert results.pvalue > 0.05
