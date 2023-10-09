"""
***Move this and copy necessary files, change remote to github***
***Create README.md***
***Create tests***

Automatically set Pdamp for Berendsen barostat by fitting to a target relaxation time, t_target.
P = P0 * exp(-t/tau) + Pset * (1 - exp(-t/tau))
P: pressure
P0: initial pressure (adjustable)
t: time
tau: time constant (adjustable)
t_set = -tau * ln(0.01) (99% of the way to set point pressure)
Pset: set point pressure
Pdamp is basically linearly related to t_target

The value of pdamp is chosen so that t_set = t_target. 
In the absence of noise or any systematic errors in fitting, pdamp is linearly related to 
t_target, so only a few short simulations are required.
"""

import argparse
import json
from pathlib import Path

import lmfit
import matplotlib.pyplot as plt
import numpy as np
from pylammpsmpi import LammpsLibrary


class SetBerendsenPdamp:
    """
    Set Pdamp for Berendsen barostat by fitting to a target relaxation time, tau.
    """

    def __init__(self, config_file):
        """
        Initialize from config_file.
        """

        with open(config_file, "r", encoding="utf-8") as fid:
            config = json.load(fid)

        # Number of cores to use for LAMMPS
        self.cores = config["CORES"]

        # Initial guess for Pdamp
        self.pdamp_initial = config["PDAMP_INITIAL"]

        # Temperature
        self.temperature = config["TSTART"]

        # Random seed
        self.seed = config["SEED"]

        # Set point pressures in LAMMPS
        self.pset = config["PSET"]

        # Simulation time for stage 2
        self.sim_time_stage2 = config["SIM_TIME_STAGE2"]

        # Target relaxation time
        self.t_target = config["T_TARGET"]

        # Tolerance for dt = |t_target - t_set|
        self.dt_tol = config["DT_TOL"]

        # Input directory
        self.indir = config["INDIR"]

        # Name of potential file to include in LAMMPS input
        self.potential_file = str(Path(self.indir, config["POTENTIAL_FILE"]))

        # LAMMPS input files (stage 1)
        if "STAGE1" in config["LAMMPS_INPUT"]:
            self.stage1_template = str(
                Path(self.indir, config["LAMMPS_INPUT"]["STAGE1"]["TEMPLATE"])
            )
            self.stage1_input = str(Path(self.indir, config["LAMMPS_INPUT"]["STAGE1"]["INPUT"]))

        # LAMMPS input files (stage 2)
        self.stage2_template = str(Path(self.indir, config["LAMMPS_INPUT"]["STAGE2"]["TEMPLATE"]))
        self.stage2_input = str(Path(self.indir, config["LAMMPS_INPUT"]["STAGE2"]["INPUT"]))

        # Output directory
        self.outdir = config["OUTDIR"]
        Path(self.outdir).mkdir(parents=True, exist_ok=True)

        # Name of data file written by LAMMPS after stage 1
        self.data_file = str(Path(self.indir, "stage1.data"))

        # Name of pressure files written by LAMMPS
        self.pressure_files = [
            str(Path(self.outdir, "pressure1.dat")),
            str(Path(self.outdir, "pressure2.dat")),
        ]

        # Time data from LAMMPS, just initialized here
        self.time = 0

        # Pressure data from LAMMPS, just initialized here
        self.pressure = 0

        # tset and pdamp data from fits
        self.pdamp = None

        # Output from lmfit
        self.fit = None

    def __call__(self):
        """
        Minimize difference between target and predicted tau by varying Pdamp.
        Plot final fit to pressure data.
        """

        self.optimize_pdamp()
        self._plot_fit()

    def optimize_pdamp(self):
        """
        Find pdamp that gives tset close to target value (dt close to 0).
        tset is nearly linearly related to pdamp, so no need to use scipy.optimize.minimize which would likely require more function evaluations (simulations).
        Stop when dt is less than self.dt_tol.
        """

        # Compute dt for self.pdamp_initial, but halve self.pdamp_initial until
        # self.sim_time_stage2 > dt + self.t_target (t_set is less than simulation time)
        # It improves convergence to start with a reasonable value of pdamp.
        dt = self.compute_dt(self.pdamp_initial)
        while dt + self.t_target > self.sim_time_stage2:
            self.pdamp_initial /= 2
            dt = self.compute_dt(self.pdamp_initial)
        self.pdamp = np.array([dt, self.pdamp_initial])
        print("dt:", dt, "| pdamp:", self.pdamp_initial)

        # Get estimate of derivative of pdamp wrt dt
        dt = self.compute_dt(1.01 * self.pdamp_initial)
        self.pdamp = np.row_stack((self.pdamp, [dt, 1.01 * self.pdamp_initial]))
        deriv = (self.pdamp[1, 1] - self.pdamp[0, 1]) / (self.pdamp[1, 0] - self.pdamp[0, 0])
        print("dt:", self.pdamp[-1, 0], "| pdamp:", self.pdamp[-1, 1])

        # Estimate pdamp that gives tset close to target
        dt1 = self.pdamp[:, 0].mean()
        pdamp1 = self.pdamp[:, 1].mean()
        pdamp2 = -deriv * dt1 + pdamp1

        # Compute dt for pdamp2 that gives dt close to 0, compute new dt
        dt = self.compute_dt(pdamp2)
        self.pdamp = np.row_stack((self.pdamp, [dt, pdamp2]))
        print("dt:", self.pdamp[-1, 0], "| pdamp:", self.pdamp[-1, 1])

        # Should be close to target, so check if difference between P0 and Pset is large enough
        # compared to fluctuations in pressure
        self._check_f()

        # Interpolate/extrapolate to get pdamp that gives dt close to 0, compute new dt
        # Fit to a second order polynomial
        # Stop when dt is less than self.dt_tol
        while abs(dt) > self.dt_tol:
            poly = np.polyfit(self.pdamp[:, 0], self.pdamp[:, 1], 2)
            pdamp = np.polyval(poly, 0)
            dt = self.compute_dt(pdamp)
            self.pdamp = np.row_stack((self.pdamp, [dt, pdamp]))
            self._check_f()
            print("dt:", self.pdamp[-1, 0], "| pdamp:", self.pdamp[-1, 1])

    def compute_dt(self, pdamp):
        """
        Compute difference between target and fit t_set.
        """

        # Edit template file
        self.edit_templates(pdamp)

        # Run LAMMPS simulations
        self.simulate()

        # Read pressure data from file
        (self.time, self.pressure) = np.loadtxt(self.pressure_files[1], unpack=True)
        self.time = self.time - self.time[0]

        # Fit to get time to set point pressure (tset) and compute dt
        dt = self._fit_tau()

        return dt

    def edit_templates(self, pdamp):
        """
        Replace [LOG_FILE], [TSTART], [SEED], [PDAMP], [PSET], [POTENTIAL_FILE], [PRESSURE_FILE],
        [SIM_TIME], & [DATA_FILE] in LAMMPS template files and write to new input files.
        """

        try:
            pdamp = pdamp[0]
        except (TypeError, IndexError):
            pass

        if self.stage1_template and self.stage1_input:
            # Stage 1 template replacements
            to_replace = [
                "[LOG_FILE]",
                "[TSTART]",
                "[SEED]",
                "[POTENTIAL_FILE]",
                "[PDAMP]",
                "[PSET]",
                "[PRESSURE_FILE]",
                "[DATA_FILE]",
            ]
            replacements = [
                str(Path(self.outdir, "stage1.log")),
                str(self.temperature),
                str(self.seed),
                self.potential_file,
                str(pdamp),
                str(self.pset[0]),
                self.pressure_files[0],
                self.data_file,
            ]

            # read stage 1 template file & make replacements
            with open(self.stage1_template, "r", encoding="utf-8") as fid:
                template_data = ''.join(fid.readlines())

            self.replace_in_template(template_data, to_replace, replacements, self.stage1_input)

        # Stage 2 template replacements
        to_replace = [
            "[LOG_FILE]",
            "[TSTART]",
            "[SEED]",
            "[DATA_FILE]",
            "[POTENTIAL_FILE]",
            "[PSET]",
            "[PDAMP]",
            "[PRESSURE_FILE]",
            "[SIM_TIME]",
        ]
        replacements = [
            str(Path(self.outdir, "stage2.log")),
            str(self.temperature),
            str(self.seed),
            self.data_file,
            self.potential_file,
            str(self.pset[1]),
            str(pdamp),
            self.pressure_files[1],
            str(self.sim_time_stage2),
        ]

        # read stage 2 template file & make replacements
        with open(self.stage2_template, "r", encoding="utf-8") as fid:
            template_data = ''.join(fid.readlines())

        self.replace_in_template(template_data, to_replace, replacements, self.stage2_input)

    def replace_in_template(self, data, toreplace, replacements, outfile):
        """
        Replace toreplace strings in data with replacements strings and write to outfile.
        """

        if len(toreplace) != len(replacements):
            raise ValueError(
                "Number of things to replace must match number of things to replace them with."
            )

        for to_rep, rep in zip(toreplace, replacements):
            data = data.replace(to_rep, rep)

        with open(outfile, "w", encoding="utf-8") as f:
            f.write(data)

    def simulate(self):
        """
        Run LAMMPS simulation(s)
        """

        # Create LammpsLibrary object
        lmp = LammpsLibrary(cores=self.cores)

        # Run stage 1 if data file doesn't exist
        if not Path(self.data_file).exists():
            str_assert = "Stage 1 input and template files must be specified in JSON input file if stage1.data file is not in INDIR."
            assert self.stage1_input and self.stage1_template, str_assert

            infile = self.stage1_input
            lmp.file(infile)
            lmp.command("clear")

            str_assert = "Your stage 1 LAMMPS input file should write a stage1.data file in INDIR."
            assert Path(self.data_file).exists(), str_assert

        # Run stage 2 simulation
        infile = self.stage2_input
        lmp.file(infile)

    def _fit_tau(self):
        """
        Fit pressure data to get tau and compute dt.
        P = P0 * exp(-t/tau) + Pset * (1 - exp(-t/tau))
        P: pressure
        P0: initial pressure (adjustable in fit)
        t: time
        tau: relaxation time (adjustable in fit)
        Pset: set point pressure
        """

        params = lmfit.Parameters()
        params.add("tau", value=self.t_target / 4, min=0.0)
        params.add("p0", value=self.pressure[0])
        params.add("pset", value=self.pset[1], vary=False)

        self.fit = lmfit.minimize(self._residual, params)

        tau = self.fit.params["tau"].value
        t_set = -tau * np.log(0.01)

        dt = t_set - self.t_target

        return dt

    def _residual(self, params):
        """
        Compute residuals for fit.
        Pactual - Ppredicted
        """

        pressure_predicted = self._pressure_function(params)
        residuals = self.pressure - pressure_predicted

        return residuals

    def _pressure_function(self, params):
        """
        (P0 * exp(-t/tau) - Pset * (1 - exp(-t/tau)))
        """

        tau = params["tau"].value
        p0 = params["p0"].value
        pset = params["pset"].value

        return p0 * np.exp(-self.time / tau) + pset * (1.0 - np.exp(-self.time / tau))

    def _check_f(self):
        """
        Check ratio of |P0 - Pset| to standard deviation of pressure, f
        """

        ind = self.time > self.t_target
        press = self.pressure[ind]
        p0 = self.fit.params["p0"].value
        dp = np.abs(p0 - self.pset[1])
        ps = press.std(ddof=1)
        f = dp / ps

        # If f is less than 2.5, raise error
        str_assert = (
            "Ratio of |P0 - Pset| to pressure fluctuations, f = "
            + str(f)
            + " Ratio of |P0 - Pset| to fluctuations in pressure is too small (< 2.5). The value of pdamp may be unreliable. Modify Pset for stage 1 or stage 2 to increase |P0 - Pset|. f should be > about 10 to get accurate values of pdamp and values < 2.5 could give unreasonable values and convergence will likely be very slow."
        )
        assert f >= 2.5, str_assert

        # If f is between 4 and 10, raise warning
        if f < 10 and f >= 2.5:
            print("\nRatio of |P0 - Pset| to pressure fluctuations, f =", f)
            print("WARNING: |P0 - Pset| is not large enough compared to fluctuations in pressure. ")
            print(
                "For f between 4 and 10, reasonable but perhaps less accurate values of pdamp will be found."
            )
            print(
                "To get more accurate values of pdamp, modify Pset for stage 1 or stage 2 to increase |P0 - Pset| above 10.\n"
            )

    def _plot_fit(self):
        """
        Save parameters to file. Plot fit to pressure data.
        """

        tau = self.fit.params["tau"].value
        tset = -tau * np.log(0.01)
        p0 = self.fit.params["p0"].value

        output = [self.temperature, self.pdamp[-1, 1], tset, tau, p0, self.pset[1]]
        output = np.array(output, dtype=float).reshape(1, -1)
        header = " ".join(["T", "pdamp", "tset", "tau", "P0", "Pset"])
        outfile = Path(self.outdir, "fit.dat")
        np.savetxt(outfile, output, header=header)

        plt.plot(self.time, self.pressure, label="data")
        plt.plot(self.time, self._pressure_function(self.fit.params), "--", label="fit")
        plt.xlabel("time")
        plt.ylabel("pressure")
        plt.legend()
        plt.title(
            "T = "
            + str(self.temperature)
            + ", Pset = "
            + str(self.pset[1])
            + ", pdamp = "
            + str(round(self.pdamp[-1, 1]))
            + ", tset = "
            + str(round(tset, 2))
            + ", tau = "
            + str(round(tau, 2))
            + ", P0 = "
            + str(round(p0))
        )
        plt.savefig(Path(self.outdir, "fit.png"))
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatically set Pdamp for Brendsen barostat by fitting to a target relaxation time, t_target."
    )
    parser.add_argument(
        "config_file",
        metavar="config_file",
        type=str,
        help="JSON file with input parameters.",
    )
    args = parser.parse_args()

    set_pdamp = SetBerendsenPdamp(args.config_file)
    set_pdamp()
