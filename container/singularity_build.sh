#!/bin/bash

# Python & LAMMPS versions from command line
# Default to 3.11 & 2Aug2023
PYTHON_VERSION=${1:-3.11}
LAMMPS_VERSION=${2:-2Aug2023}

# Output directory
outdir=/home/bnovak1/singularity/set_berendsen_pdamp

# Build singularity image from docker image
singularity build ${outdir}/python_${PYTHON_VERSION}_lammps_${LAMMPS_VERSION,,}.sif \
                  docker-daemon:python_${PYTHON_VERSION}_lammps_${LAMMPS_VERSION,,}:latest
ln -s -n ${outdir}/python_${PYTHON_VERSION}_lammps_${LAMMPS_VERSION,,}.sif