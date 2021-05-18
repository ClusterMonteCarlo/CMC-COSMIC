.. _examples:

##################
Initial Conditions
##################

Running CMC comes in two distinct stages: generating the initial conditions, 
and then running CMC on those initial conditions.  We briefly cover here how to 
create initial conditions and how to run and restart CMC. 

=============================
Point-mass Initial Conditions
=============================

.. ipython::
   :suppress:

   In [143]: import os

   In [144]: from matplotlib.pylab import *

   In [145]: style.use(os.path.abspath('./source/matplotlib_style_file'))

   In [146]: ion()

Generating initial conditions is the first step in running CMC.  Since Version 
3.4, COSMIC has been able to create CMC initial conditions, using a combination 
of the native binary samplers and dynamical samplers for drawing positions and 
velocities of clusters in virial and hydrostatic equilibrium.  COSMIC can 
currently generate clusters from Plummer, Elson, or King profiles, as detailed below.  

The examples below are also found in the examples folder of the main GitHub 
repository.  Note that there are two distinct options for CMC clusters: 
``cmc_point_mass`` (which we use in the Plummer and Elson examples below) will 
produce a cluster of point-mass particles: no stellar radii, mass functions, or 
binaries, typical for purely dynamics problems.  You can see it's options with:

.. ipython:: python

        from cosmic.sample.sampler import cmc

        help(cmc.get_cmc_point_mass_sampler) 


The more general ``cmc`` sampler (which we use in the King profile example 
below) has all the additional features for generating stellar binaries that 
COSMIC has.  You can see it's (significantly more diverse) options with

.. ipython:: python

        from cosmic.sample.sampler import cmc

        help(cmc.get_cmc_sampler) 

Plummer Sphere
--------------

One of the most common initial conditions for a star cluster are those of 
`Plummer (1911) 
<https://ui.adsabs.harvard.edu/abs/1911MNRAS..71..460P/abstract>`_.  They are 
not the most representitive of star clusters in the local universe (see the 
King profile below), their popularity endures because their main quantities 
(such as enclosed mass, density, velocity dispersion, etc.) can be expressed analytically.

To generate a Plummer profile with 10,000 particles (here we only consider 
point-mass particles, as opposed to real stars), we can use the following 
COSMIC functions.  Note that this returns two HDF5 tables, containing all the 
relevant parameters for the single stars and binaries.  The two tables are 
linked, such 
that the ``Singles`` table contains all the single stars and binary centers of mass, as well as indicies pointing to the relevant binaries in ``Binaries``.  

.. ipython:: python
    :okwarning:

        from cosmic.sample import InitialCMCTable

        Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='plummer', size=10000, r_max=100)

        print(Singles) 

We then scale the cluster to `Hénon units 
<https://ui.adsabs.harvard.edu/abs/2014arXiv1411.4936H/abstract>`_, where 
:math:`G = 1`, the cluster mass is :math:`M_{c}=1`, the cluster energy is 
:math:`E=-0.25`, and the virial radius is :math:`r_v \equiv G M_c^2 / 2|U| = 
1`.  Note that is automatically done in ``InitialCMCTable.write`` if the 
cluster isn't already normalized, but we do it here explicitly.

.. ipython:: python

        InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries) 

Finally, we can save the initial conditions to an HDF5 file (for use in CMC) with:

.. ipython:: python

        InitialCMCTable.write(Singles, Binaries, filename="plummer.hdf5") 


We can check that the Plummer function reproduces what we would expect from 
analytic predictions.  The enclosed mass a plummer sphere is given by

.. math::

   M(r) = M_{\rm total}\left(1 + \frac{a^2}{r^2}\right)^{-3/2}

where :math:`a` is an arbitrary scale factor (which is :math:`3\pi / 16` when 
the virial radius is normalized to 1).  If we compare the mass-weighted 
cumulative radii of our ``Singles`` Pandas table to the analytic results, we 
can see:

.. ipython:: python

        import numpy as np

        import matplotlib.pyplot as plt

        r_grid = np.logspace(-1.5,2,100)

        m_enc = (1 + 1/r_grid**2)**-1.5 

        rv = 16/(3*np.pi) # virial radius for a Plummer sphere

        plt.plot(r_grid/rv,m_enc,lw=3);

        plt.hist(Singles.r,weights=Singles.m,cumulative=True,bins=r_grid);

        plt.xscale('log')

        plt.xlabel("Radii",fontsize=15);

        plt.ylabel(r"$M (< r) / M_{\rm total}$",fontsize=15);

        @savefig plot_plummer.png width=6in
        plt.legend(("Plummer","COSMIC Samples"),fontsize=14);



Elson Profile
--------------
The `Elson (1987) 
<https://ui.adsabs.harvard.edu/abs/1987ApJ...323...54E/abstract>`_ profile is a 
generalization of the Plummer profile that has been shown to better represent 
young massive clusters in the local universe.  The density at a 
radius :math:`r` is given by 

.. math::

   \rho(r) = \rho_{0}\left(1 + \frac{a^2}{r^2}\right)^{-\frac{\gamma + 1}{2}}

Note that :math:`\gamma = 4` gives a Plummer profile (the above code actually 
just calls the Elson profile generator with :math:`\gamma=4`), though most 
young clusters are better fit with :math:`\gamma\sim2-3`.  The enclosed mass is 
correspondingly more complicated:

.. math::

   M(r) = \frac{4 \pi \rho_0}{3} r^3 \,_2F_1\left(\frac{3}{2},\frac{\gamma + 1}{2} ; \frac{5}{2} ; -\frac{r^2}{a^2}\right)  

Where :math:`\,_2F_1` is the ordinary hypergeometric function.  

Unlike both the Plummer and King profiles, the distribution function for the 
Elson profile cannot be written analytically.  To genereate the initial 
conditions, we directly integrate the density and potential functions to 
numerically compute :math:`f(E)`, and draw our velocity samples from that (see 
appendix B of `Grudic et al., 2018 
<https://ui.adsabs.harvard.edu/abs/2018MNRAS.481..688G/abstract>`_).  This 
produces a handful 
of warnings in the SciPy integrators, but the profiles that it generates are correct.

To generate an Elson profile with :math:`\gamma=3`, we can use

.. ipython:: python
    :okwarning:

    Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='elson', gamma=2.5, size=10000, r_max=100) 


Comparing with the theoretical calculation for the enclosed mass, we find similarly good agreement:

.. ipython:: python

        from scipy.special import hyp2f1

        gamma = 2.5

        def m_enc(gamma,r,rho_0):
            return 4*np.pi*rho_0/3 * r**3 * hyp2f1(1.5,(gamma+1)/2,2.5,-r**2)

        rv = 6. ## virial radius for a gamma=2.5 Elson profile

        rho_0 = 1/m_enc(gamma,100*rv,1)

        plt.plot(r_grid/rv,m_enc(gamma,r_grid,rho_0),lw=3); # note we scale by rv, rather than set scale factor

        plt.hist(Singles.r,weights=Singles.m,cumulative=True,bins=r_grid);

        plt.xscale('log')

        plt.xlabel("Radii",fontsize=15);

        plt.ylabel(r"$M (< r) / M_{\rm total}$",fontsize=15);

        @savefig plot_elson.png width=6in
        plt.legend(("Elson ($\gamma=2.5$)","COSMIC Samples"),fontsize=14);

King Profile
------------
An idealized cluster in thermodynamic equilibrium could be described as an 
isothermal sphere, where the velocities of stars resembled a Maxwell-Boltzmann 
distribution.  But the isothermal sphere has infinite mass, and in any 
realistic star cluster, the distribution of stars should go to zero near the 
tidal boundary.  The `King (1966) <https://ui.adsabs.harvard.edu/abs/1966AJ.....71...64K/abstract>`_ profile 
acomplishes this by sampling from a lowered isothermal distribution

.. math::

        f(E) = \begin{cases}
            \rho_0 (2\pi\sigma^2)^{-3/2}(e^{E/\sigma^2})&    E > 0;\\
            0& E \leq 0
        \end{cases}

The King initial conditions can be created with COSMIC using:

.. ipython:: python
    :okwarning:

    Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='king', w_0=6, size=10000, r_max=100)

The analytic form of :math:`M(<r)` cannot be written down for a King profile, 
but we can solve the ODE directly (this is done when generating the samples)

.. ipython:: python

        from cosmic.sample.cmc import king
        radii,rho,phi,m_enc = king.integrate_king_profile(6)

        rho /= m_enc[-1]
        m_enc /= m_enc[-1]

        rv = king.virial_radius_numerical(radii, rho, m_enc) # just compute R_v numerically
        plt.plot(radii/rv,m_enc,lw=3);
        plt.hist(Singles.r,weights=Singles.m,cumulative=True,bins=r_grid);
        plt.xscale('log')
        plt.xlim(0.05,10);
        plt.xlabel("Radii",fontsize=15);
        plt.ylabel(r"$M (< r) / M_{\rm total}$",fontsize=15);

        @savefig plot_king.png width=6in
        plt.legend(("King ($w_0=6$)","COSMIC Samples"),fontsize=14);



.. _realistic-initial-conditions:

============================
Realistic Initial Conditions
============================


So far the above examples have only used the ``cmc_point_mass`` sampler.  To 
generate `realistic` initial conditions, with stellar masses and binaries, we 
want to use the ``cmc`` sampler instead.  This enables all the additional 
options found in the ``independent`` population sampler that COSMIC uses `(see 
here for more details) 
<https://cosmic-popsynth.github.io/COSMIC/runpop/index.html#independent>`_.  

To generate the above King profile, but with all the additional stellar physics (initial mass function, binaries, etc.) we would use

.. ipython:: python

        Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.1, primary_model='kroupa01',
                                                    ecc_model='thermal', porb_model='log_uniform', qmin=-1.0,
                                                    cluster_profile='king', met=0.00017, size=100000,w_0=6,
                                                    params=os.path.abspath('../../examples/KingProfile.ini'),
                                                    seed=12345,virial_radius=1,tidal_radius=1e6)

This example is also found in the `examples 
<https://github.com/ClusterMonteCarlo/CMC-COSMIC/tree/master/examples>`_ folder 
in the CMC repository.  Note that unlike the above examples, here we require an 
ini file (or a BSEDict), to supply the necessary stellar physics options for 
generating stars and binaries with 
BSE.
