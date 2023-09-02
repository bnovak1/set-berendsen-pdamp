import matplotlib.pyplot as plt
import numpy as np

infile = "Pdamp_vs_P0.dat"
data = np.loadtxt(infile)
data = data[data[:, 0].argsort()]

outfile = "Pdamp_vs_P0.png"
plt.plot(data[:, 0]/6398.903051160765, data[:, 1], "o-")
plt.xlabel("$f$")
plt.ylabel("$pdamp$")
plt.savefig(outfile, bbox_inches="tight")