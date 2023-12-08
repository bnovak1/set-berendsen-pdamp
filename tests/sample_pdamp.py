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
    fit_file = Path(base_dir, "output/fit.dat")

    data_output = pd.DataFrame([], columns=["seed", "pdamp", "t_set", "tau", "P0"])

    for seed in seeds:
        with open(infile, "r", encoding="utf-8") as jf:
            config = json.load(jf)

        config["SEED"] = int(seed)

        with open(infile, "w", encoding="utf-8") as jf:
            json.dump(config, jf, indent=4)

        set_pdamp = SetBerendsenPdamp(infile)
        set_pdamp()

        fit_df = pd.read_csv(fit_file, sep="\s+", header=0)
        cols = fit_df.columns
        fit_df.drop(columns=[cols[-1]], inplace=True)
        fit_df.columns = cols[1:]
        fit_df = pd.concat(
            [pd.DataFrame([seed], columns=["seed"]), fit_df[["pdamp", "t_set", "tau", "P0"]]],
            axis=1,
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
