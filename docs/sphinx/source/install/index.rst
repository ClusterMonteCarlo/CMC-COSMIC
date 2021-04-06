.. _install:

############
Installation
############

=======================
Generic HPC Instuctions
=======================
To build CMC on a standard high-performance computing (HPC) system, you'll need a few things beforehand.  These should be standard on just about any cluster environment; if they're not, try using the Docker image that is linked to below!

The specific requirements are

 * Fortran and C compilers
 * Message Passing Interface (MPI)
 * Gnu Scientific Library (GSL)
 * HDF5 Libraries
 * cmake (version 3.12 or higher)
 * cfitsio (optional; only needed for legacy initial condition generators)

These should be available on any HPC system.  If they're available, we suggest using the Intel version of the compilers and MPI.  On the machines we've tested, this produces a significant speedup over GCC/OpenMPI. 

Once you have either installed or module loaded the above requirements, then download the CMC-COSMIC package with:

.. code-block:: bash 
   
    git clone git@github.com:ClusterMonteCarlo/CMC-COSMIC.git --recurse-submodules

GCC
_____
Using the MPI-wrapped commands for your C and Fortran compilers, build the into the bin folder with

.. code-block:: bash

    cd CMC-COSMIC
    mkdir build bin
    cd build
    FC=mpifort CC=mpicc cmake .. -DCMAKE_INSTALL_PREFIX=../bin 
    make install

Intel
_____
If you are using the intel compilers, you can instead use

.. code-block:: bash

    cd CMC-COSMIC
    mkdir build bin
    cd build
    FC=mpif77 CC=mpiicc cmake .. -DCMAKE_INSTALL_PREFIX=../bin 
    make install

====================
Specific HPC Systems 
====================

Quest
-----
FOR SCOTTY

Bridges
-------
On Bridges 2 (the XSEDE machine) start by importing cmake and the intel compilers

.. code-block:: bash

   module load intelmpi/20.4-intel20.4

Then follow the **intel** instructions above 

Vera 
----
On Vera (McWilliams center machine) or Henon (Rodriguez group machine), first load

.. code-block:: bash

    module load intel/20.2
    module load cmake/3.18.1

Then follow the **intel** instructions above 

================
Docker Container
================
FOR SCOTTY (+ Carl)
