.. _install:

############
Installation
############

=======================
Generic HPC Instuctions
=======================
To build CMC on a standard high-performance computing (HPC) system, you'll need 
a few things beforehand.  These should be standard on just about any cluster 
environment; if they're not, try using the Docker image that is linked to 
below!

The specific requirements are

 * Fortran and C compilers
 * Message Passing Interface (MPI)
 * Gnu Scientific Library (GSL)
 * HDF5 Libraries
 * cmake (version 3.12 or higher)
 * cfitsio (optional; only needed for legacy initial condition generators)

These should be available on any HPC system.  If they're available, we suggest 
using the Intel version of the compilers and MPI.  On the machines we've 
tested, this produces a significant speedup over GCC/OpenMPI. 

Once you have either installed or module loaded the above requirements, then 
download the CMC-COSMIC package with:

.. code-block:: bash 
   
    git clone https://github.com/ClusterMonteCarlo/CMC-COSMIC.git --recurse-submodules

GCC
_____
Using the MPI-wrapped commands for your C and Fortran compilers, build the into the bin folder with

.. code-block:: bash

    cd CMC-COSMIC
    mkdir build
    cd build
    FC=mpifort CC=mpicc cmake .. -DCMAKE_INSTALL_PREFIX=../CMC 
    make install

.. note::

    By default, the above command will put all the CMC executables, libraries, and include files into a folder called `CMC` in the main repository folder.  You 
    can customize where this is installed by changing the ``-DCMAKE_INSTALL_PREFIX`` option above

Intel
_____
If you are using the intel compilers, you can instead use

.. code-block:: bash

    cd CMC-COSMIC
    mkdir build
    cd build
    FC=mpiifort CC=mpiicc cmake .. -DCMAKE_INSTALL_PREFIX=../CMC 
    make install

=================
Installing COSMIC
=================

There are several ways to install COSMIC along with CMC.  If you have a version 
of ipython installed (preferably though anaconda) that you can run ``pip 
install`` with, then cmake can install the version of COSMIC that comes with 
CMC (and all the associated dependencies).  Simply add the 
``-DBUILD_COSMIC=ON`` flag to the cmake step in the installation:

.. code-block:: bash

    FC=mpiifort CC=mpiicc cmake .. -DBUILD_COSMIC=ON -DCMAKE_INSTALL_PREFIX=../CMC

This will install the version of COSMIC that is included with CMC (in the ./src folder) into your python path.

You can also use the stand alone version of COSMIC that is available `here 
<https://cosmic-popsynth.github.io/COSMIC/install/index.html>`_

====================
Specific HPC Systems 
====================
Specific instructions for several of the HPC systems used by the CMC/COSMIC 
teams are below (as well as some XSEDE machines).  If you are interested in 
adding more instructions, just let us know!

Quest
_____
To build Quest with Intel

.. code-block:: bash

    module load cmake/3.15.4 hdf5/1.10.7-openmpi-4.0.5-intel-19.0.5.281 mpi/openmpi-4.0.5-intel-19.0.5.281 gsl/2.5-intel-19.0.5.281

Then follow the **intel** instructions above

Bridges
_______
On Bridges 2 (the XSEDE machine) start by importing cmake and the intel compilers

.. code-block:: bash

   module load intelmpi/20.4-intel20.4

Then follow the **intel** instructions above 

Vera/Hénon
__________
On Vera (McWilliams center machine) or Hénon (Rodriguez group machine), first load

.. code-block:: bash

    module load intel/20.2
    module load cmake/3.18.1

Then follow the **intel** instructions above 

Frontera 
________
On Frontera (the new fancy TACC machine from XSEDE) the commands:

.. code-block:: bash

    module load TACC intel impi hdf5 gsl
    
then following the **intel** instructions above should work (h/t to Mike Grudić for testing CMC there)

==========
Containers
==========

Singularity
-----------
To use a container to run CMC on Quest, you must use singularity. First you must pull the container from docker locally

.. code-block:: bash
    singularity pull docker://clustermontecarlo/cmc

In order to run a container which has a MPI-enabled program on an HPC system, you must use an OpenMPI program on the host system whose major version (think 4.X.X) matches that used inside the container. In this case, we compile CMC against OpenMPI 4.X.X. For those running on a system using SLURM, a submission script might look at follows

.. code-block:: bash
    #!/bin/bash 
    #SBATCH -J container_test
    #SBATCH --error=error.out
    #SBATCH --output=output.out
    #SBATCH --nodes=2
    #SBATCH --ntasks-per-node=4
    #SBATCH --mem=0
    #SBATCH --time=00:15:00
    #SBATCH --account=XXXXX
    #SBATCH --partition=XXXX 

    module purge all
    module load mpi/openmpi-4.0.5-gcc-10.2.0
    module load singularity

    mpirun -np ${SLURM_NTASKS} singularity exec -B /projects:/projects openmpi.sing /opt/local/bin/cmc CMC_Params.ini initial

Docker
------
.. code-block:: bash

    docker pull clustermontecarlo/cmc
