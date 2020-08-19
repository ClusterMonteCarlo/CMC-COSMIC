# Installation Instructions

```
mkdir build
cd build
CONDA_BUILD_SYSROOT=/ FC=mpifort CC=mpicc cmake .. -DCMAKE_INSTALL_PREFIX=.
CONDA_BUILD_SYSROOT=/ make install
```
