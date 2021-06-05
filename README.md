# The Cluster Monte Carlo Code, CMC-COSMIC
A Henon style orbit-averaging code for collisional stellar dynamics, with all the necessary physics for modeling the long-term evolution of dense star clusters.  

See the documentation [here](https://clustermontecarlo.github.io/)

# Installation Instructions (See [here](https://clustermontecarlo.github.io/CMC-COSMIC/install/index.html) for more details)

Checkout the code with (be sure to add the `--recurse-submodules` flag to download COSMIC as well!)

```
git clone https://github.com/ClusterMonteCarlo/CMC-COSMIC.git --recurse-submodules
```

You'll need MPI, HDF5, GSL, and cmake (version 3.12 or higher) to compile the code.  These should be available on any HPC system (or easily installed with your favorite package manager on Linux).  Once downloaded, compile the code with

```
cd CMC-COSMIC
mkdir build
cd build
FC=mpifort CC=mpicc cmake .. -DCMAKE_INSTALL_PREFIX=../CMC
make install
```

This will create a build and bin folder in the top level of the CMC-COSMIC package with the relevant executables.

# Citing the code
If you use CMC-COSMIC in your research, please cite [Joshi et al., 2000](https://ui.adsabs.harvard.edu/abs/2000ApJ...540..969J/abstract), [Pattabiraman et al., 2013](https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract), and [Rodriguez et al., in prep](https://i.pinimg.com/originals/ca/fe/32/cafe32d3992642c2f21b56d87d8b4768.gif), as well as the COSMIC paper [Breivik et al., 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract).  In a pinch, you can just cite Rodriguez et al., but please cite all 4 if possible.

# Badges
[![Unit Test CMC](https://github.com/ClusterMonteCarlo/CMC-COSMIC/actions/workflows/python-package.yml/badge.svg)](https://github.com/ClusterMonteCarlo/CMC-COSMIC/actions/workflows/python-package.yml)
[![Docker](https://img.shields.io/docker/v/clustermontecarlo/cmc.svg)](https://hub.docker.com/repository/docker/clustermontecarlo/cmc)

# Use CMC in a container

## Docker
```
docker pull clustermontecarlo/cmc
```
## Singularity
```
singularity pull docker://clustermontecarlo/cmc
```
