# set-berendsen-pdamp

The response in the pressure due to the [Berendsen barostat](http://www.sklogwiki.org/SklogWiki/index.php/Berendsen_barostat) is related to the bulk modulus ([compressibility](http://www.sklogwiki.org/SklogWiki/index.php/Compressibility#Isothermal_compressibility)). For the [press/berendsen fix](https://docs.lammps.org/fix_press_berendsen.html) in [LAMMPS](https://www.lammps.org/#gsc.tab=0), you can either compute the bulk modulus of your system and use the modulus keyword and pdamp or use the default value of modulus and set pdamp to get the desired response.

This package automatically sets pdamp for the Berendsen barostat in LAMMPS by fitting to a target time to the set point pressure, $t_{target}$. For a given value of pdamp, the pressure versus time data is fit to:

$$ P = P_0 \exp(-t/\tau) + P_{set} (1 - \exp(-t/\tau)) $$

$P$: pressure  
$t$: time  
$P_{set}$: set point pressure  
$P_0$: initial pressure (adjustable in fit)  
$\tau$: time constant (adjustable in fit)  
$t_{set} = -\tau \ln(0.01)$ (99% of the way from $P_0$ to $P_{set}$)

The value of pdamp is chosen so that $t_{set} \approx t_{target}$. In the absence of noise or any systematic errors in fitting, pdamp is linearly related to $t_{target}$, so only a few short simulations are required.

## Installation

### This package

### LAMMPS & pylammpsmpi

If you don't need any extra LAMMPS packages and you are using conda, then installing [pylammpsmpi](https://pylammpsmpi.readthedocs.io/en/latest/installation.html) will also install LAMMPS.

```shell
conda install pylammpsmpi -c conda-forge
```

Installing pylammpsmpi using pip will not install LAMMPS. If you are using pip for all of your python packages, then you must compile LAMMPS yourself even if you don't need any extra LAMMPS packages.

```shell
pip install pylammpsmpi
```

If you need any extra LAMMPS packages, you must compile LAMMPS yourself and use pip to install pylammpsmpi. LAMMPS must be installed as a library and with the python module. The [calphy package documentation](https://calphy.org/en/latest/gettingstarted.html) explains how to do this. Also, see the [LAMMPS documentation for installing with the python module](https://docs.lammps.org/Python_install.html). Note that you should install any other LAMMPS packages you need using additional `-D` flags with the `cmake` command.

## Usage

### Inputs

#### JSON input file

#### LAMMPS input file(s), stage1.data

#### Choosing set point pressure(s)

![](README/Pdamp_vs_f.png)

$f=\frac{\left|P_0 - P_{set}\right|}{\sigma_P}$

### Python script

```python
from set_Berendsen_Pdamp import SetBerendsenPdamp
setter = SetBerendsenPdamp('./config.json')
setter()
```

### Command line

```shell
python set_Berendsen_Pdamp.py ./config.json
```
