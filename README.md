# Installation Instructions

## Create virtual environment
###Linux
```
conda install -c conda-forge openmpi-mpicc gcc_osx-64 gxx_linux-64 gsl cfitsio gfortran_linux-64 python=3.7 autoconf libtool automake pkg-config binutils hdf5
```

###MacOS
```
conda install -c conda-forge openmpi-mpicc clang_osx-64 clangxx_osx-64 gsl cfitsio gfortran_osx-64 python=3.7 autoconf libtool automake pkg-config binutils hdf5
```

## install CMC
```
mkdir build
cd build
CONDA_BUILD_SYSROOT=/ FC=mpifort CC=mpicc cmake .. -DCMAKE_INSTALL_PREFIX=.
CONDA_BUILD_SYSROOT=/ make install
```
