.. _examples:

############
Getting Started
############

Running CMC comes in two distinct stages: generating the initial conditions, and then running CMC on those initial conditions.  We briefly cover here how to 
create initial conditions and how to run and restart CMC. 

==================
Initial Conditions
==================

.. ipython::
   :suppress:

   In [144]: from matplotlib.pylab import *

   In [145]: style.use('./matplotlib_style_file')

   In [146]: ion()

Generating initial conditions is the first step in running CMC.  Since Version 3.4, COSMIC has been able to create CMC initial conditions, using a combination 
of the native binary samplers and dynamical samplers for drawing positions and velocities of clusters in virial and hydrostatic equilibrium.  COSMIC can 
currently generate clusters from Plummer, Elson, or King profiles, as detailed below.  

The examples below are also found in the examples folder of the main GitHub repository.  Note that there are two distinct options for CMC clusters: 
``cmc_point_mass`` (which we use in the Plummer and Elson examples below) will produce a cluster of point-mass particles: no stellar radii, mass functions, or binaries, typical for purely dynamics problems.  You can see it's options with:

.. ipython::

    In [2]: from cosmic.sample.sampler import cmc

    In [3]: help(cmc.get_cmc_point_mass_sampler)


The more general ``cmc`` sampler (which we use in the King profile example below) has all the additional features for generating stellar binaries that COSMIC 
has.  You can see it's (significantly more diverse) options with

.. ipython::

    In [2]: from cosmic.sample.sampler import cmc

    In [3]: help(cmc.get_cmc_sampler)

Plummer Sphere
--------------

One of the most common initial conditions for a star cluster are those of `Plummer (1911) <https://ui.adsabs.harvard.edu/abs/1911MNRAS..71..460P/abstract>`_.  
They are not the most representitive of star clusters in the local universe (see the King profile below), their popularity endures because their main quantities 
(such as enclosed mass, density, velocity dispersion, etc.) can be expressed analytically.

To generate a Plummer profile with 10,000 particles (here we only consider point-mass particles, as opposed to real stars), we can use the following COSMIC 
functions.  Note that this returns two HDF5 tables, containing all the relevant parameters for the single stars and binaries.  The two tables are linked, such 
that the ``Singles`` table contains all the single stars and binary centers of mass, as well as indicies pointing to the relevant binaries in ``Binaries``.  

.. ipython::
    :okwarning:

    In [1]: from cosmic.sample import InitialCMCTable
    
    In [2]: Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='plummer', size=10000, r_max=100)

    In [3]: print(Singles)

We then scale the cluster to `Henon units <https://ui.adsabs.harvard.edu/abs/2014arXiv1411.4936H/abstract>`_, where :math:`G = 1`, the cluster mass is :math:`M_{c}=1`, the cluster energy is :math:`E=-0.25`, and the virial radius is :math:`r_v \equiv G M_c^2 / 2|U| = 1`.  Note that is automatically done in ``InitialCMCTable.write`` if the cluster isn't already normalized, but we do it here explicitly.

.. ipython::

    In [3]: InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries)

Finally, we can save the initial conditions to an HDF5 file (for use in CMC) with:
.. ipython::

    In [4]: InitialCMCTable.write(Singles, Binaries, filename="plummer.hdf5")


We can check that the Plummer function reproduces what we would expect from analytic predictions.  The enclosed mass a plummer sphere is given by

.. math::

   M(r) = M_{\rm total}\left(1 + \frac{a^2}{r^2}\right)^{-3/2}

where :math:`a` is an arbitrary scale factor (which we set to :math:`3\pi / 16` when the virial radius is normalized to 1).  If we compare the mass-weighted 
cumulative radii of our ``Singles`` Pandas table to the analytic results, we can see:

.. ipython::

    In [3]: import numpy as np

    In [4]: import matplotlib.pyplot as plt

    In [5]: a = 3*np.pi/16

    In [6]: r_grid = np.logspace(-1.5,1.5,100)

    In [7]: m_enc = (1 + a**2/r_grid**2)**-1.5 

    In [8]: plt.plot(r_grid,m_enc,lw=3);

    In [9]: plt.hist(Singles.r,weights=Singles.m,cumulative=True,bins=r_grid);

    In [10]: plt.xscale('log')
    
    In [11]: plt.xlabel("Radii",fontsize=15);

    In [12]: plt.ylabel(r"$M (< r) / M_{\rm total}$",fontsize=15);

    @savefig plot_simple.png width=7in
    In [13]: plt.legend(("Theory","COSMIC Samples"),fontsize=14);



Elson Profile
------------
The `Elson (1987) <https://ui.adsabs.harvard.edu/abs/1987ApJ...323...54E/abstract>`_ profile is a generalization of the Plummer profile that has been shown to 
better represent young massive clusters in the local universe.  The density at a 
radius :math:`r` is given by 

.. math::

   \rho(r) = \rho_{0}\left(1 + \frac{a^2}{r^2}\right)^{-\frac{\gamma + 1}{2}}

Note that :math:`\gamma = 4` gives a Plummer profile (the above code actually just calls the Elson profile generator with :math:`\gamma=4`), though most young 
clusters are better fit with :math:`\gamma\sim2-3`.  The enclosed mass is correspondingly more complicated:

.. math::

   M(<r) = \frac{4 \pi \rho_0}{3} r^3 \,_2F_1\left(\frac{3}{2},\frac{\gamma + 1}{2} ; \frac{5}{2} ; -\frac{r^2}{a^2}\right)  

Where :math:`\,_2F_1` is the ordinary hypergeometric function.  

Unlike both the Plummer and King profiles, the distribution function for the Elson profile cannot be written analytically.  To genereate the initial conditions, 
we directly integrate the density and potential functions to numerically compute :math:`f(E)`, and draw our velocity samples from that (`see appendix B of Grudic et al., 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.481..688G/abstract>`_).  This produces a handful 
of warnings in the SciPy integrators, but the profiles that it generates are correct.

To generate an Elson profile with :math:`\gamma=3`, we can use
.. python::

    from cosmic.sample import InitialCMCTable
    
    Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='elson', gamma=3, size=10000, r_max=100)

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
