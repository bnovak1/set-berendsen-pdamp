"""
Automatically set Pdamp for Berendsen barostat by fitting to a target relaxation time, `t_target`.
`P = P0 * exp(-t/tau) + Pset * (1 - exp(-t/tau))`
`P`: pressure
`P0`: initial pressure (adjustable)
`t`: time
`tau`: time constant (adjustable)
`t_set = -tau * ln(0.01)` (99% of the way to set point pressure)
`Pset`: set point pressure
`Pdamp` is basically linearly related to `t_target`

The value of `Pdamp` is chosen so that `t_set = t_target`. 
In the absence of noise or any systematic errors in fitting, `Pdamp` is linearly related to 
`t_target`, so only a few short simulations are required.

The class `SetBerendsenPdamp` is used to optimize the value of `Pdamp` for the Berendsen barostat. The class takes a JSON configuration file as input.

1. Initialize the class with the input parameters from the JSON file.
2. Optimize the value of `Pdamp` by minimizing the difference between the target relaxation time and the predicted relaxation time. Each of the following steps is repeated for each step of the optimization until the difference between the target and fit relaxation times is below a specified tolerance:
    1. Edit LAMMPS input template files to replace placeholders with the appropriate values of `Pdamp` and other parameters.
    2. Run LAMMPS simulations using the edited input files.
    3. Read the pressure data from the simulations and fit it to obtain the relaxation time.
    4. Compare the fitted relaxation time to the target relaxation time and adjust the value of `Pdamp` accordingly.
3. Plot the fit to the pressure data.

The value of `Pdamp` is stored in the `pdamp` attribute of the class. The fit parameters, including the relaxation time, are saved to a file. A plot of the fit to the pressure data is also saved as an image file.

Example Usage:
config_file = "input.json"
set_pdamp = SetBerendsenPdamp(config_file)
set_pdamp()
"""

import argparse
from collections.abc import Iterable
import json
from pathlib import Path
import warnings

import lmfit
import numpy as np
from pylammpsmpi import LammpsLibrary


class SetBerendsenPdamp:
    """
    Set Pdamp for Berendsen barostat by fitting to a target relaxation time, t_target.

    Attributes:
    -----------
    cores : int
        Number of cores to use for LAMMPS simulations
    pdamp_initial : float
        Initial guess for Pdamp
    temperature : float
        Temperature
    seed : int
        Random seed
    pset : list
        Set point pressures in LAMMPS
    sim_time_stage2 : float
        Simulation time for stage 2
    t_target : float
        Target relaxation time
    dt_tol : float
        Tolerance for dt = |t_target - t_set|
    indir : str
        Input directory
    potential_file : str
        Name of potential file to include in LAMMPS input
    stage1_template : str
        LAMMPS input file name (stage 1) template
    stage1_input : str
        LAMMPS input file name (stage 1)
    stage2_template : str
        LAMMPS input file name (stage 2) template
    stage2_input : str
        LAMMPS input file name (stage 2)
    outdir : str
        Output directory name
    data_file : str
        Name of data file written by LAMMPS after stage 1
    pressure_files : list
        Name of pressure files written by LAMMPS
    time : float
        Time data read from LAMMPS pressure_files. Initialized to 0.
    pressure : float
        Pressure data read from LAMMPS pressure_files. Initialized to 0.
    pdamp : float
        tset and pdamp data from fits. Initialized to None.
    fit : lmfit object
        Output from lmfit. Initialized to None.
    """

    def __init__(self, config_file):
        """
        Initialize the SetBerendsenPdamp object from a JSON configuration file.

        Args:
        -----------
        config_file : str
            Name of the JSON configuration file.
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

        # Make sure stage1.data file is in INDIR if stage 1 files are not specified
        if not Path(self.data_file).exists():
            if not hasattr(self, "stage1_template") or not hasattr(self, "stage1_input"):
                raise AttributeError(
                    "Stage 1 input and template files must be specified in the JSON input file if the stage1.data file is not in INDIR."
                )

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
        Minimize the difference between target and predicted relaxation time by varying Pdamp.
        Save the final fit to pressure data.
        """
        self.optimize_pdamp()
        self._save_fit()

    def optimize_pdamp(self):
        """
        Find pdamp that gives tset close to the target value (dt close to 0).
        tset is nearly linearly related to pdamp, so there is no need to use scipy.optimize.minimize which would likely require more function evaluations (simulations).
        Stop when dt is less than self.dt_tol.
        """

        # Compute dt for self.pdamp_initial, but halves self.pdamp_initial until
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
        Edit the LAMMPS template files, run the LAMMPS simulation(s), read the pressure data, and fit the data to get the difference between the target and fitted values of the relaxation time, dt.

        Args:
        -----------
        pdamp : float
            The value of pdamp to use in the LAMMPS simulation(s).

        Returns:
        -----------
        dt : float
            The computed difference between the target and fitted values of the relaxation time.
        """

        # Edit template file
        self.edit_templates(pdamp)

        # Run LAMMPS simulation(s)
        if not Path(self.data_file).exists():
            self.simulate(1)
        self.simulate(2)

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

        Args:
        -----------
        pdamp : float
            Value of Pdamp to use in LAMMPS input files.
        """

        if isinstance(pdamp, Iterable):
            if len(pdamp) > 0:
                pdamp = pdamp[0]
            else:
                raise ValueError("pdamp must be a scalar or an iterable with at least one value.")

        # If stage 1 input and template files are specified, make replacements
        if self.stage1_template and self.stage1_input:
            # Stage 1 template replacements
            replacements = {
                "[LOG_FILE]": str(Path(self.outdir, "stage1.log")),
                "[TSTART]": str(self.temperature),
                "[SEED]": str(self.seed),
                "[POTENTIAL_FILE]": self.potential_file,
                "[PDAMP]": str(pdamp),
                "[PSET]": str(self.pset[0]),
                "[PRESSURE_FILE]": self.pressure_files[0],
                "[DATA_FILE]": self.data_file,
            }

            # read stage 1 template file & make replacements
            with open(self.stage1_template, "r", encoding="utf-8") as fid:
                template_data = "".join(fid.readlines())

            self.replace_in_template(template_data, replacements, self.stage1_input)

        # Stage 2 template replacements
        replacements = {
            "[LOG_FILE]": str(Path(self.outdir, "stage2.log")),
            "[TSTART]": str(self.temperature),
            "[SEED]": str(self.seed),
            "[DATA_FILE]": self.data_file,
            "[POTENTIAL_FILE]": self.potential_file,
            "[PSET]": str(self.pset[1]),
            "[PDAMP]": str(pdamp),
            "[PRESSURE_FILE]": self.pressure_files[1],
            "[SIM_TIME]": str(self.sim_time_stage2),
        }

        # read stage 2 template file & make replacements
        with open(self.stage2_template, "r", encoding="utf-8") as fid:
            template_data = "".join(fid.readlines())

        self.replace_in_template(template_data, replacements, self.stage2_input)

    def replace_in_template(self, data, replacements, outfile):
        """
        Replace keys in data with values from replacements and write to outfile.

        Args:
        -----------
        data : str
            The string to perform replacements on.
        replacements : dict
            The keys are strings to replace in data with values from replacements.
        outfile : str
            Name of LAMMPS file to write to after replacements are made.
        """
        for to_rep, rep in replacements.items():
            data = data.replace(to_rep, rep)

        with open(outfile, "w", encoding="utf-8") as f:
            f.write(data)

    def simulate(self, stage_number):
        """
        Run LAMMPS simulation(s)

        Args:
        -----------
        stage_number: int
            The stage number of the simulation to run. Must be 1 or 2.

        Raises:
        -----------
        ValueError
            If stage_number is not 1 or 2.
        """

        # Create LammpsLibrary object
        lmp = LammpsLibrary(cores=self.cores)

        # Run simulation
        if stage_number == 1:
            infile = self.stage1_input
            lmp.file(infile)

            str_assert = "Your stage 1 LAMMPS input file should write a stage1.data file in INDIR."
            if not Path(self.data_file).exists():
                raise FileNotFoundError(str_assert)

        elif stage_number == 2:
            infile = self.stage2_input
            lmp.file(infile)

        else:
            raise ValueError("Stage_number must be 1 or 2.")

    def _fit_tau(self):
        """
        Fit pressure data to get tau and compute dt.

        This function fits the pressure vs. time data to obtain the time constant (tau) which is used to compute the relaxation time (t_set) and the difference between t_set and the target relaxation time (t_target):

        P = P0 * exp(-t/tau) + Pset * (1 - exp(-t/tau))

        where:
        - P: pressure
        - P0: initial pressure (adjustable in fit)
        - t: time
        - tau: time constant (adjustable in fit)
        - Pset: set point pressure (not adjustable in fit)

        t_set = -tau * ln(0.01) (99% of the way to set point pressure)

        Returns:
        -----------
        dt : float
            Absolute difference between t_set and its target value, t_target.
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
        Compute residuals for fit of pressure (P) vs. time (t) data to P = P0 * exp(-t/tau) + Pset * (1 - exp(-t/tau))

        Args:
        -----------
        params : lmfit.Parameters
            lmfit parameters (tau, P0, Pset) for the pressure function.

        Returns:
        --------
        residuals : numpy.ndarray
            Array of residuals between actual and predicted pressure values.
        """
        pressure_predicted = self._pressure_function(params)
        residuals = self.pressure - pressure_predicted

        return residuals

    def _pressure_function(self, params):
        """
        Fitted pressure as a function of time using P = P0 * exp(-t/tau) + Pset * (1 - exp(-t/tau))

        Args:
        -----------
        params : lmfit.Parameters
            lmfit parameters (tau, P0, Pset) for the pressure function.

        Returns:
        -----------
        pfit : numpy.ndarray
            The fitted pressure as a function of time.
        """

        tau = params["tau"].value
        p0 = params["p0"].value
        pset = params["pset"].value

        pfit = p0 * np.exp(-self.time / tau) + pset * (1.0 - np.exp(-self.time / tau))

        return pfit

    def _check_f(self):
        """
        Check ratio of |P0 - Pset| to standard deviation of pressure, f

        This function calculates the ratio of the absolute difference between the initial pressure and the target pressure to the standard deviation of the pressure. If the ratio is less than 2.5, an error is raised, as the value of pdamp may be unreliable. If the ratio is between 4 and 10, a warning is printed, as the value of pdamp may be less accurate.
        To get more accurate values of pdamp, modify Pset for stage 1 or stage 2 or increase the system size (reduces standard deviation) to increase the ratio above 10.
        """

        ind = self.time > self.t_target
        press = self.pressure[ind]
        p0 = self.fit.params["p0"].value
        dp = np.abs(p0 - self.pset[1])
        ps = press.std(ddof=1)
        f = dp / ps

        # If f is less than 2.5, raise error
        str_error = (
            f"The ratio of |P0 - Pset| to pressure fluctuations is f = {f}. The ratio of |P0 - Pset| to fluctuations in pressure is too small (< 2.5). The value of pdamp may be unreliable. Modify Pset for stage 1 or stage 2 to increase |P0 - Pset|. f should be > about 10 to get accurate values of pdamp and values < 2.5 could give unreasonable values and convergence will likely be very slow."
        )
        if f < 2.5:
            raise ValueError(str_error)

        # If f is between 4 and 10, raise warning
        str_warn = (
            f"The ratio of |P0 - Pset| to pressure fluctuations is f = {f}. For f between 4 and 10, reasonable but perhaps less accurate values of pdamp will be found. To get more accurate values of pdamp, modify Pset for stage 1 or stage 2 to increase |P0 - Pset| above 10."
        )
        if f < 10 and f >= 2.5:
            warnings.warn(UserWarning(str_warn))
        

    def _save_fit(self):
        """
        Saves temperature, Pdamp, t_set, tau, P0, & Pset; time, pressure, and the fit to the pressure data.
        """
        tau = self.fit.params["tau"].value
        t_set = -tau * np.log(0.01)
        p0 = self.fit.params["p0"].value

        output = {}
        output["temperature"] = self.temperature
        output["pdamp"] = self.pdamp[-1, 1]
        output["t_set"] = t_set
        output["tau"] = tau
        output["P0"] = p0
        output["Pset"] = self.pset[1]
        output["time"] = self.time
        output["pressure"] = self.pressure
        output["fit"] = self._pressure_function(self.fit.params)
        
        with open(Path(self.outdir, "fit.json"), "w", encoding="utf-8") as jf:
            json.dump(output, jf, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatically set Pdamp for Brendsen barostat by fitting to a target relaxation time, t_target."
    )
    parser.add_argument("config_file", type=str, help="Path to JSON file with input parameters.")
    args = parser.parse_args()

    set_pdamp = SetBerendsenPdamp(args.config_file)
    set_pdamp()
