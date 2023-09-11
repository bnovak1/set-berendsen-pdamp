"""
Plot Pdamp vs |P0 - Pset|/fluctuations to show that Pdamp is less accurate when P0 
is closer to the size of the pressure fluctuations.
"""

import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Standard deviation of pressure data from some simulations with f > 3
dir = str(Path("pressure_data/pressure*.dat"))
press_files = glob.glob(dir)

for press_file in press_files:
    
    p = np.loadtxt(press_file, skiprows=1)
    ind = p[:, 0] >= 1.0 # Time > t_target
    p = p[ind, :]
    
    try:
        press = np.row_stack((press, p))
    except NameError:
        press = p

ps = press.std(ddof=1)

infile = "Pdamp_vs_P0.dat"
data = np.loadtxt(infile)
data = data[data[:, 0].argsort()]
f = data[:, 0]/ps

ind = f > 8.5
plt.plot(f[ind], data[ind, 1], "go", mfc="none")
ind = f <= 8.5
plt.plot(f[ind], data[ind, 1], "ro", mfc="none")
ylims = plt.ylim()
plt.plot([8.5, 8.5], ylims, "k--")
plt.ylim(ylims)
plt.xlim(2.5, plt.xlim()[1])

plt.xlabel("$f$")
plt.ylabel("$pdamp$")

outfile = "Pdamp_vs_f.png"
plt.savefig(outfile, bbox_inches="tight")
plt.close()