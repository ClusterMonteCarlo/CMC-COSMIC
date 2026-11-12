.. _parallel:

############
Code Structure and Parallelization
############

There are various routines which have varying dependencies and accesses between 
elements of the data structures. To minmize inter-process communication 
required by these routines as much as possible, we partition the data such that 
the number of stars held by each processor is a multiple of ``MIN_CHUNK_SIZE`` 
(input parameter, defaults to 20, refer to `Pattabiraman et al. (2013) 
<https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_  for a 
justification for the choice of this value, and detailed explanation). The 
remainder of stars (which is less than ``MIN_CHUNK_SIZE``) after division goes 
to the last processor.


Most routines do computations on a star by star basis, and access only the 
properties of the star being processed. However, the values of gravitational 
potential, masses of the stars and radial positions of the stars are accessed 
in a complex, data-dependent manner throughout the code. Hence, we store these 
values in separate arrays and duplicate/copy these across all processors. Since 
these properties change during a timestep, we synchronize these arrays at the 
end of each time step.

Like mentioned above, most routines perform computations on a star by star 
basis. The only routine that differs from this pattern is the sorting routine 
that sorts all the stars by their radial positions. We use a parallel sorting 
algorithm called Sample Sort. 

The parallel sorting routine does the sorting and automatically redistributes 
the data to the processors, however, there is no control over the number of 
stars that will end up on each processor. So, after the sorting takes place, we 
exchange data between processors to adhere to the data partitioning scheme 
mentioned above.

Programming Guidelines for Developers
-------------------------------------
It is very necessary to read `Pattabiraman et al. (2013) 
<https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_ before 
proceeding. The paper covers most of the code and parallelization aspects which 
in itself serves as a guide for development. Here, we cover some more low level 
coding details to aid a developer who wants to modify the code. An introduction 
to MPI might also be very useful, and in fact necessary for any non-trivial 
code development. 

The parallelization follows a SIMD model. Any line of code is executed by all 
processors. Unlike a sequential code, variables might have different values on 
different processors. A simple example is the variable ``myid`` which stores 
the ID of the processor. If p processors are used, the IDs go from ``0`` to 
``p-1``, and each processor will contain the value of its ID in this variable. 
On the other hand the variable ``procs`` stores the total number of processors 
used. One can see that the variable ``procs`` will have the same value on all 
processors whereas ``myid`` will have different value depending on the 
processor.

The data/stars are initially partitioned, after reading from the input file, 
into the user specified number of processors. The partitioning scheme is 
described in detail in the paper, and is chosen in such a way as to minimize 
communication required during the program execution. Binaries are partitioned 
among processors as well, and are indexed similar to the serial code i.e. the 
variable ``binind`` of a star, if greater than 0, refers to the index of the 
corresponding binary in the binary array.

.. _bharat: http://adsabs.harvard.edu/abs/2013ApJS..204...15P
__ bharat_
The variables ``mpiBegin``, ``mpiEnd``, and the arrays ``Start`` and ``End`` 
are used to implement and store the data partitioning scheme. These are used in 
almost every routine of the code. Let us try to understand these with an 
example. Say, we start with an initial of 10000 stars, and 4 processors. The 
data is partitioned, and each processor has 2500 stars. These are stored in the 
star array. Please note that in each timestep all the stars are sorted by their 
radial position. So, positions of all stars in processor 0 are lesser than 
those in processor 1 which in turn are lesser than the ones in processor 2 and 
so on. Now, each local array is indexed from 1 to 2500. However, if you imagine 
an array containing the entire set of 10000 stars, the indices of the stars 
local to a processor will have a different index, which will depend on the 
processor ID, in this imaginary global array. For instance the local stars 1 to 
2500 in processor 1 will have a "global" index 2501 to 5000. This "global" 
index is require at various points in the code since some routines need some 
quantities of the entire set of stars. These quantities (mass, positions and 
potential) are duplicated across all processors, please see `Pattabiraman et 
al. (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_ 
for a detailed description. The variables ``mpiBegin`` and ``mpiEnd`` store the 
global index of the first and the last star in a given processor. For our 
example, in processor 0, these will be 1 and 2500 (the 0th star is always left 
out as it is marked as a sentinel); for processor 1 these will be 2501 and 5000 
and so on. The ``Start`` and ``End`` arrays store these values of all the 
processors. ``Start[myid]`` essentially is the same as ``mpiBegin``, and 
``End[myid]`` has the same value as ``mpiEnd``. In addition there are the 
``mpiDisp`` and ``mpiLen`` arrays which also store essentially the same 
information but as displacement offsets and lengths, but these are used only at 
very few places.

The total number of stars in a given processor can be obtained by 
``mpiEnd-mpiBegin+1``. ``clus.N_MAX`` represents the current total number of 
stars in the simulation, whereas ``clus.N_MAX_NEW`` is the total number of 
local stars in the processor i.e. ``mpiEnd-mpiBegin+1``. The function 
``mpiFindIndicesCustom()`` sets ``mpiBegin`` and ``mpiEnd``, whereas 
``findLimits()`` populates the ``Start`` and ``End`` arrays. When one iterates 
over the local number of stars, and given a star ``i``, one can get its global 
index using the function ``get_global_idx(i)``.

The random number generation is parallelized, and all details are hidden, so a 
developer does not have to mess with the details as long as one follows a 
similar rng calls as an already existing one.

Example
-----
While introducing any new code, think in parallel. For instance say you are 
calculating some cumulative quantity like:


.. code-block:: c 

    double Q
    for(i=1; i<=clus.N_MAX_NEW; i++) {
        Q = star[i].a + star[i].b;
    }


This might appear as if it's correct, although when this runs in parallel it is 
far from it. ``Q`` will only contain the local sum ``Q`` for each processor. 
After this one needs to aggregate ``Q`` across all processors. For this you'll 
have to use either ``MPI_Reduce`` or ``MPI_Allreduce`` depending on whether you 
want the aggregated quantity ``Q`` on one or all of the processors. The above 
piece of code can be fixed by:

.. code-block:: c 

   double Q, local_Q;

   for(i=1; i<=clus.N_MAX_NEW; i++) {
       local_Q = star[i].a + star[i].b;
   }

   MPI_Allreduce(&local_Q, &Q, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

Here we store the local value of ``Q`` into a temporary variable ``local_Q``, 
after which we aggregate the partial sums across all processors using the 
``MPI_Allreduce`` call which makes sure all processors receive the aggregated 
value in ``Q``. Here, we use ``MPI_COMM_WORLD``, which is the default 
communicator that includes all processors.

This is probably the simplest use of parallel programming and use of MPI. For 
more complex operations which involve complex data dependencies such as 
communicating neighbors/ghost particles, gathering, scattering data etc., one 
might need to have more expertise in parallel programming using MPI. So, a 
tutorial on MPI is strongly advised before adding any non-trivial piece of 
code. The most common MPI calls used in the code are 
``MPI_Reduce/MPI_Allreduce``, ``MPI_Gather/MPI_Allgather``, ``MPI_Alltoall, 
MPI_Bcast``, ``MPI_Scan/MPI_Exscan``. The signatures of all MPI calls can be 
found here: url{http://www.mpich.org/static/docs/v3.1/www3/}

I/O
-----
Output is another part which one might need to add, and is non-trivial to do in 
parallel. However, we have introduced a decently good framework to do parallel 
output so the programmer doesn't have to delve into the MPI IO details. If a 
file needs to be written by only one node, it's simple. Just do:

.. code-block:: c 

   if(myid == <NODE_ID>) {
      FILE *fp;
      fopen(...);
      fprintf(...);
      fclose(...);
   }

Sometimes, the data that needs to be written out might be distributed across 
nodes. Provided these are ASCII data files with file sizes of a few KBs to MBs 
(such as log files, and not snapshots which write out nearly entire data) one 
can use MPI IO to write them out in parallel. An example is as follows:

.. code-block:: c 

   MPI_File mpi_binfp;
   char mpi_binfp_buf[10000], mpi_binfp_wrbuf[10000000];
   long long mpi_binfp_len=0, mpi_binfp_ofst_total=0;
   sprintf(filename, "a_e2.%04ld.dat", tcount);
   MPI_File_open(MPI_COMM_MC, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_binfp);
   MPI_File_set_size(mpi_binfp, 0);

   for (j=1; j<=clus.N_MAX; j++) {
      if (star_array[j].binind) {
         parafprintf(binfp, "%g %g\n", binary_array[star_array[j].binind].a, sqr(binary_array[star_array[j].binind].e));
      }
   }

   mpi_para_file_write(mpi_binfp_wrbuf, &mpi_binfp_len, &mpi_binfp_ofst_total, &mpi_binfp);
   MPI_File_close(&mpi_binfp);

The files are opened by all processors using MPI-IO. Each processor writes the 
data into a string/char buffer. At the end of the timestep, all processors 
flush the data from the buffers into the corresponding files in parallel using 
MPI-IO. The code uses 5 variables for this process - the MPI-IO file pointer, 
which follows the format ``mpi_<fptr_name>``, 2 char buffers, which have the 
format ``mpi_<fptr_name>_buf`` and ``mpi_<fptr_name>_wrbuf``, and an 
int/longlong variables to maintain the length of the buffer (format 
``mpi_<fptr_name>_len``) and the offset in the file (format 
``mpi_<fptr_name>_ofst_total``) where data has to be written. Given all these 
variables, a call to the function ``parafprintf(<fptr_name>, <args>)`` which 
follows a pattern similar to ``fprintf``, but underneath writes the given args 
to the respective buffers and updates the offsets etc. Once all writing has 
been done, each processor has different data in their respective buffers which 
needs to be flushed out into the file. A call to ``mpi_para_file_write()`` 
takes care of this.


