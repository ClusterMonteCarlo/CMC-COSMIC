.. _examples:

############
Examples
############

FOR CARL

==================
Initial Conditions
==================

Generating initial conditions is the first step in running CMC.  Since Version 3.4, COSMIC has been able to create CMC initial conditions, using a combination 
of the native binary samplers and dynamical samplers for drawing positions and velocities of clusters in virial and hydrostatic equilibrium.  COSMIC can 
currently generate clusters from Plummer, Elson, or King profiles, as detailed below.  

To see the arguments necessary to call the CMC sampler use the help function:

.. ipython::
   :suppress:

   In [144]: from matplotlib.pylab import *

   In [145]: ion()

.. ipython::

    In [2]: from cosmic.sample.sampler import cmc

    In [3]: help(cmc.get_cmc_sampler)

The examples below are also found in the examples folder of the main GitHub repository.

Plummer Sphere
--------------

One of the most common initial conditions for a star cluster are those of `Plummer (1911) <https://ui.adsabs.harvard.edu/abs/1911MNRAS..71..460P/abstract>`_.  
They are not the most representitive of star clusters in the local universe (see the King profile below), their popularity endures because their main quantities 
(such as enclosed mass, density, velocity dispersion, etc.) can be expressed analytically.

To generate a Plummer profile with 10,000 particles (here we only consider point-mass particles, as opposed to real stars), we can use the following COSMIC 
functions.  Note that this returns two HDF5 tables, containing all the relevant parameters for the single stars and binaries.  The two tables are linked, such 
that the ``Singles`` table contains all the single stars and binary centers of mass, as well as indicies pointing to the relevant binaries in ``Binaries``.  

.. ipython::

    In [1]: from cosmic.sample import InitialCMCTable
    
    In [2]: Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='plummer', size=10000, r_max=100)

    In [3]: print(Singles)

We then scale the cluster to `Henon units <https://ui.adsabs.harvard.edu/abs/2014arXiv1411.4936H/abstract>`_, where :math:`G = 1`, the cluster mass is :math:`M_{c}=1`, the cluster energy is :math:`E=-0.25`, and the virial radius is :math:`r_v \equiv G M_c^2 / 2|U| = 1`.  Note that is automatically done in ``InitialCMCTable.write`` if the cluster isn't already normalized, but we do it here explicitly.

.. code::

    In [3]: InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries)

Finally, we can save the initial conditions to an HDF5 file (for use in CMC) with:
.. code::

    In [4]: InitialCMCTable.write(Singles, Binaries, filename="plummer.hdf5")


We can check that the Plummer function reproduces what we would expect from analytic predictions.  The enclosed mass a plummer sphere is given by

.. math::

   M(r) = M_{\rm total}\left(1 + \frac{a^2}{r^2}\right)^{-3/2}

where :math:`a` is an arbitrary scale factor (which we set to :math:`3\pi / 16` when the virial radius is normalized to 1).  If we compare the mass-weighted 
cumulative radii of our ``Singles`` Pandas table to the analytic results, we can see:

.. ipython::

    In [3]: import numpy as np
    In [3]: import matplotlib.pyplot as plt

    In [3]: a = 3*np.pi/16
    In [3]: r_grid = np.logspace(-1.5,1.5,100)
    In [3]: m_enc = (1 + a**2/r**2)**-1.5 
    In [3]: plt.plot(r_grid,M_enclosed,lw=2)

    In [3]: plt.hist(Singles.r,weights=Singles.m,cumulative=True,bins=r_grid);
    In [3]: plt.xscale('log')

    @savefig plot_simple.png width=5in
    In [3]: plt.legend(("Theory","COSMIC Samples"))



Elson Profile
------------
The `Elson (1987) <https://ui.adsabs.harvard.edu/abs/1987ApJ...323...54E/abstract>`_ profile is a generalization of the Plummer profile.  

King Profile
------------
same shit here...

===========
Running CMC
===========
put some example C code here


==============
Restarting CMC
==============

`
