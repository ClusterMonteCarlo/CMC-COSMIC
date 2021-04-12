.. _intro:

############
Introduction
############

Code History
=======================================
This code offers one approach for solving the gravitational `N`-body problem. The basic problem is as follows: you have a system of `N` particles (stars, in our case) they are all interacting with one another via the gravitational force. How does the system evolve over time? There are several ways to solve the `N`-body problem, each method having its pros and cons. The two most common methods in use today are the directy `N`-body method and the Monte Carlo method. CMC uses the Monte Carlo method of `Henon (1971) <https://ui.adsabs.harvard.edu/abs/1971Ap%26SS..14..151H/abstract>`_. We will discuss them both briefly before moving into a more detailed discussion of our MC technique.

In the direct `N`-body method, the equations of motion for each star are solved directly by summing the gravitational force of every star on every star. Each star feels a gravitational pull from all N-1 other stars, and in each timestep, their positions and velocities are updated according to the total acceleration force it feels over that timestep. This is the most exact method of solving the `N`-body problem, but it has its costs: the computation time scales as `N`:sup:`2`. The Monte Carlo method represents an approximate way to solve the `N`-body problem that makes use of three assumptions: first, that the dynamical evolution is dominated by two-body relaxation, and second, that the cluster has a sufficient number of particles that relaxation time of most stars is significantly longer than their orbital timescales.  

These two assumptions are often referred to as the Fokker-Planck approximation, and are appropriate for large clusters like the globular clusters observed in our own galaxy (and the younger super star clusters observed in nearby galaxies).  It is this limit that the Henon approach (and its descendants, such as CMC) operates.  Our implementation further assumeps spherical symmetry and orbit averaging (like Henon's original implementation), which allows us to model old spherical star clsuters significantly faster than direct `N`-body integrators. Unlike Henon's original implementation, modern computer hardware and distributed memory parallelization allows CMC to create star-by-star realizations of the cluster, where each dynamical particle is an individual star, with as much additional physics (e.g. stellar evolution, strong dynamical encounters, binary formation, etc.) as the user desires. Integration is done on a relaxation timescale, with the code scaling as O(N Log N) in the number of particles.

CMC has been developed over many years starting with `Joshi et al. (2000 <https://ui.adsabs.harvard.edu/abs/2000ApJ...540..969J/abstract>`_, `2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...550..691J/abstract>`_, `Fregeau et al. (2003) <https://ui.adsabs.harvard.edu/abs/2003ApJ...593..772F/abstract>`_.  It has been upgraded to include all the relevant physics for modeling dense spherical star clusters, such as strong dynamical encounters (`Fregeau & Rasio 2007 <https://ui.adsabs.harvard.edu/abs/2007ApJ...658.1047F/abstract>`_), single and binary stellar evoltuion (`Chatterjee et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010ApJ...719..915C/abstract>`_), central massive black holes (`Umbreit et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...750...31U/abstract>`_), three-body binary formation (`Morscher et al., 2014 <https://ui.adsabs.harvard.edu/abs/2013ApJ...763L..15M/abstract>`_), relativistic dynamics (`Rodriguez et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018PhRvL.120o1101R/abstract>`_) and more.  CMC has been parallelized using the Message Passing Interface (MPI), allowing for it's use on both large workstations and distributed-memory high-performance computing clusters (`Pattabiraman et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_), and has been shown to produce similar results to state-of-the-art direct `N`-body simulations (`Rodriguez et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2109R/abstract>`_).

The original implementation of stellar evolution in CMC used the original binary stellar evolution package of `Hurley et al. (2002) <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_, with updates from `Chatterjee et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010ApJ...719..915C/abstract>`_ and `Rodriguez et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018PhRvL.120o1101R/abstract>`_.  The public release of CMC presented here is pinned to the `COSMIC <https://cosmic-popsynth.github.io/>`_ package for binary population synthesis (`Breivik et al. 2019 <https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract>`_).  This symbiosis is performed at a repository level: any changes to COSMIC are automatically propogated to CMC.  


Flowchart
=========================

Here we show the basic flow of decisions made during every CMC timestep.  Each of these steps are described in one (or several of the CMC papers below).  As an 
initial first guess, most of the relevant physics is described and properly cited in `Pattabiraman et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_. Otherwise, the refences below are   

1. `Joshi et al. (2000 <https://ui.adsabs.harvard.edu/abs/2000ApJ...540..969J/abstract>`_, `2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...550..691J/abstract>`_
2. `Fregeau & Rasio 2007 <https://ui.adsabs.harvard.edu/abs/2007ApJ...658.1047F/abstract>`_ 
3. `Chatterjee et al. 2010 <https:/ui.adsabs.harvard.edu/abs/2010ApJ...719..915C/abstract>`_
4. `Goswami et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...752...43G/abstract>`_ 
5. `Pattabiraman et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_
6. `Breivik et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract>`_
7. `Morscher et al., 2014 <https://ui.adsabs.harvard.edu/abs/2013ApJ...763L..15M/abstract>`_

.. figure:: even_bigger_flowchart.jpg

   Flowchart of the main dynamics loop in CMC 
