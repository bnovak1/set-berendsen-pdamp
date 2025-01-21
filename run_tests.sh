#!/bin/bash

# Run tests in singularity containers for different Python & LAMMPS versions

# Flag to run test_optimization or not. Default to 1 (run)
opt_flag=${1:-1}

python_version="3.11"
lammps_version="2Aug2023"

# Run tests not requiring LAMMPS    
singularity exec container/python_${python_version}_lammps_${lammps_version,,}.sif \
    pytest -k "not test_optimization"

# Run test_optimization which requires LAMMPS if opt_flag=1
# For some reason this test never finishes when run with the other tests
# Use --capture=no to see output
if [ "$opt_flag" -eq 1 ]; then
    singularity exec container/python_${python_version}_lammps_${lammps_version,,}.sif \
        pytest --capture=no -k test_optimization
fi
    
