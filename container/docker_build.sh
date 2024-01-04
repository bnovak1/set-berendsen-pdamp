#!/bin/bash

# Python & LAMMPS version from command line, default to 3.11 & 2Aug2023
PYTHON_VERSION=${1:-3.11}
LAMMPS_VERSION=${2:-2Aug2023}

echo "Building docker image python_${PYTHON_VERSION}_lammps_${LAMMPS_VERSION,,}"
docker build -t python_${PYTHON_VERSION}_lammps_${LAMMPS_VERSION,,} \
             --build-arg PYTHON_VERSION=${PYTHON_VERSION} \
             --build-arg LAMMPS_VERSION=${LAMMPS_VERSION} \
             -f Dockerfile .



