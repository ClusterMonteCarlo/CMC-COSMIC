.. _intro:

############
Introduction
############

Code Philosophy 
=======================================
This code offers one approach for solving the gravitational `N`-body problem. 
The basic problem is as follows: you have a system of `N` particles (stars, in 
our case) they are all interacting with one another via the gravitational 
force. How does the system evolve over time? There are several ways to solve 
the `N`-body problem, each method having its pros and cons. The two most common 
methods in use today are the directy `N`-body method and the Monte Carlo 
method. CMC uses the Monte Carlo method of `Henon (1971) 
<https://ui.adsabs.harvard.edu/abs/1971Ap%26SS..14..151H/abstract>`_. We will 
discuss them both briefly before moving into a more detailed discussion of our 
MC technique.

In the direct `N`-body method, the equations of motion for each star are solved 
directly by summing the gravitational force of every star on every star. Each 
star feels a gravitational pull from all N-1 other stars, and in each timestep, 
their positions and velocities are updated according to the total acceleration 
force it feels over that timestep. This is the most exact method of solving the 
`N`-body problem, but it has its costs: the computation time scales as 
`N`:sup:`2`. The Monte Carlo method represents an approximate way to solve the 
`N`-body problem that makes use of three assumptions: first, that the dynamical 
evolution is dominated by two-body relaxation, and second, that the cluster has 
a sufficient number of particles that relaxation time of most stars is 
significantly longer than their orbital timescales.  

These two assumptions are often referred to as the Fokker-Planck approximation, 
and are appropriate for large clusters like the globular clusters observed in 
our own galaxy (and the younger super star clusters observed in nearby 
galaxies).  It is this limit that the Henon approach (and its descendants, such 
as CMC) operates.  Our implementation further assumeps spherical symmetry and 
orbit averaging (like Henon's original implementation), which allows us to 
model old spherical star clsuters significantly faster than direct `N`-body 
integrators. Unlike Henon's original implementation, modern computer hardware 
and distributed memory parallelization allows CMC to create star-by-star 
realizations of the cluster, where each dynamical particle is an individual 
star, with as much additional physics (e.g. stellar evolution, strong dynamical 
encounters, binary formation, etc.) as the user desires. Integration is done on 
a relaxation timescale, with the code scaling as :math:`O(N \log N)` in the 
number of particles.


Flowchart
=========================

Here we show the basic flow of decisions made during every CMC timestep.  Each 
of these steps are described in one (or several of the CMC papers below).  As 
an initial first guess, most of the relevant physics is described and properly 
cited in `Pattabiraman et al. 2013 
<https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_. Otherwise, 
the refences below are   

1. `Joshi et al. (2000 <https://ui.adsabs.harvard.edu/abs/2000ApJ...540..969J/abstract>`_, `2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...550..691J/abstract>`_
2. `Fregeau & Rasio 2007 <https://ui.adsabs.harvard.edu/abs/2007ApJ...658.1047F/abstract>`_ 
3. `Chatterjee et al. 2010 <https:/ui.adsabs.harvard.edu/abs/2010ApJ...719..915C/abstract>`_
4. `Goswami et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...752...43G/abstract>`_ 
5. `Pattabiraman et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_
6. `Breivik et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract>`_
7. `Morscher et al., 2014 <https://ui.adsabs.harvard.edu/abs/2013ApJ...763L..15M/abstract>`_

.. figure:: even_bigger_flowchart.jpg

   Flowchart of the main dynamics loop in CMC 
