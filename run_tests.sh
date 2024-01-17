#!/bin/bash

# Run tests in singularity containers for different Python & LAMMPS versions

# Flag to run test_optimization or not. Default to 0 (don't run)
opt_flag=${1:-0}

# Use all Python & LAMMPS versions if opt_flag=1
# Otherwise, use only one LAMMPS version (2Aug2023)
if [ "$opt_flag" -eq 1 ]; then
    python_versions=("3.8" "3.9" "3.10" "3.11")
    lammps_versions=("23Jun2022" "2Aug2023")
else
    python_versions=("3.8" "3.9" "3.10" "3.11")
    lammps_versions=("2Aug2023")
fi

# Loop over Python & LAMMPS versions
for python_version in "${python_versions[@]}"; do
    for lammps_version in "${lammps_versions[@]}"; do

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
    
    done
done