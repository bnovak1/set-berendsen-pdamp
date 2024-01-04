# Start with a base ubuntu image
FROM ubuntu:22.04

# Set the working directory
WORKDIR /app

# Update and upgrade
RUN apt-get update && apt-get upgrade -y

# Install necessary ubuntu packages
RUN apt-get install -y cmake curl gcc g++ openssh-client

# Install miniconda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /app/miniconda3 \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && /app/miniconda3/bin/conda clean --all --yes

# Add miniconda/bin to path
ENV PATH=/app/miniconda3/bin:$PATH

# Activate conda environment
RUN conda init \
    && bash

# Add conda-forge to channels
RUN conda config --append channels conda-forge

# Install python packages
ARG PYTHON_VERSION=3.11
RUN conda install -y python=$PYTHON_VERSION \
pytest lmfit mpi4py numpy openmpi pandas scipy \
    && conda clean --all --yes

# Copy LAMMPS source code
# Pre-downloaded with curl -LO https://download.lammps.org/tars/lammps-${LAMMPS_VERSION}.tar.gz
ARG LAMMPS_VERSION=2Aug2023
COPY lammps-${LAMMPS_VERSION}.tar.gz /app
RUN tar -xzf lammps-${LAMMPS_VERSION}.tar.gz \
    && rm lammps-${LAMMPS_VERSION}.tar.gz \
    && mv lammps-${LAMMPS_VERSION} lammps

# Install LAMMPS as a library and python module
RUN cd lammps \
    && mkdir build_lib \
    && cd build_lib \
    && cmake -D BUILD_LIB=ON -D BUILD_SHARED_LIBS=ON -D BUILD_MPI=ON -D PKG_MPIIO=ON -D LAMMPS_EXCEPTIONS=yes -D PKG_MANYBODY=ON -D PKG_MISC=ON -D PKG_EXTRA-COMPUTE=ON -D PKG_EXTRA-DUMP=ON -D PKG_EXTRA-FIX=ON -D PKG_EXTRA-PAIR=ON ../cmake \
    && make -j 4 \
    && cp liblammps${SHLIB_EXT}* ../src \
    && cd ../src \
    && make install-python \
    && mkdir -p $CONDA_PREFIX/include/lammps \
    && cp library.h $CONDA_PREFIX/include/lammps \
    && cp liblammps${SHLIB_EXT}* $CONDA_PREFIX/lib

# Install pylammpsmpi
RUN pip install pylammpsmpi

# Set the user & run conda init to add to .bashrc of user
RUN groupadd -r user && useradd -r -g user -G sudo -m -s /bin/bash user
USER user
WORKDIR /home/user
RUN conda init