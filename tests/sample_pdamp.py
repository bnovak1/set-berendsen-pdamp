"""
Run multiple optimizations with different seeds to get a distribution of pdamp values to use for testing.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from set_berendsen_pdamp import SetBerendsenPdamp


def multi_sims(seeds, base_dir="tests", infile="input/config.json"):
    """
    Run multiple optimizations with different seeds to get a distribution of pdamp values
    """

    infile = Path(base_dir, infile)
    fit_file = Path(base_dir, "output/fit.json")

    with open(infile, "r", encoding="utf-8") as jf:
        config = json.load(jf)

    data_output = pd.DataFrame([], columns=["seed", "pdamp", "t_set", "tau", "P0"])

    for seed in seeds:
        config["SEED"] = int(seed)

        with open(infile, "w", encoding="utf-8") as jf:
            json.dump(config, jf, indent=4)

        set_pdamp = SetBerendsenPdamp(infile)
        set_pdamp()

        with open(fit_file, "r", encoding="utf-8") as jf:
            fit_output = json.load(jf)

        fit_df = pd.DataFrame(
            {
                "seed": [seed],
                "pdamp": [fit_output["pdamp"]],
                "t_set": [fit_output["t_set"]],
                "tau": [fit_output["tau"]],
                "P0": [fit_output["P0"]],
            }
        )

        data_output = pd.concat([data_output, fit_df], axis=0)

    return config, data_output


if __name__ == "__main__":
    config, data_output = multi_sims(np.unique(np.random.randint(1e6, 1e7, size=1000)))

    output = {
        "description": "Values of pdamp obtained with different random seeds. The config key is the input file used for the optimization except the seed is changed. The data key contains the seed, pdamp, t_set, tau, and P0 values obtained from the optimizations.",
        "config": config,
        "data": data_output.to_dict(orient="list"),
    }
    with open(Path("tests/pdamp_samples.json"), "w", encoding="utf-8") as jf:
        json.dump(output, jf, indent=4)
