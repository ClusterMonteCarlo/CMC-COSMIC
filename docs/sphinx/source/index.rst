.. CMC documentation master file, created by
   sphinx-quickstart on Mon Aug  3 21:43:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================================
The Cluster Monte Carlo Code, CMC
=================================
Welcome to the documentation for CMC, an `N`-body code for collisional stellar dynamics based on the original method of `Henon (1971) <https://ui.adsabs.harvard.edu/abs/1971Ap%26SS..14..151H/abstract>`_.  CMC has been under development for more than 20 years, starting with `Joshi et al. (2000 <https://ui.adsabs.harvard.edu/abs/2000ApJ...540..969J/abstract>`_, `2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...550..691J/abstract>`_, `Fregeau et al. (2003) <https://ui.adsabs.harvard.edu/abs/2003ApJ...593..772F/abstract>`_.  It has been upgraded to include all the relevant physics for modeling dense spherical star clusters, such as strong dynamical encounters (`Fregeau & Rasio 2007 <https://ui.adsabs.harvard.edu/abs/2007ApJ...658.1047F/abstract>`_), single and binary stellar evoltuion (`Chatterjee et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010ApJ...719..915C/abstract>`_), central massive black holes (`Umbreit et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...750...31U/abstract>`_), three-body binary formation (`Morscher et al., 2014 <https://ui.adsabs.harvard.edu/abs/2013ApJ...763L..15M/abstract>`_), relativistic dynamics (`Rodriguez et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018PhRvL.120o1101R/abstract>`_) and more.  CMC is parallelized using the Message Passing Interface (MPI), allowing it to run on workstations and distributed-memory high-performance computing clusters (`Pattabiraman et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJS..204...15P/abstract>`_). It has been shown to produce similar results to state-of-the-art direct `N`-body simulations (`Rodriguez et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2109R/abstract>`_).

The public version of CMC is pinned to the `COSMIC <https://cosmic-popsynth.github.io/>`_ package for binary population synthesis (`Breivik et al. 2019 <https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract>`_), which itself was originally based on the version of BSE (`Hurley et al. 2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_) that was developed in CMC over many years (`Chatterjee et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010ApJ...719..915C/abstract>`_, `Rodriguez et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016PhRvD..93h4029R/abstract>`_, and `2018 <https://ui.adsabs.harvard.edu/abs/2018PhRvL.120o1101R/abstract>`_).  Both projects have advanced significantly since then, but are pinned at the repository level: COSMIC is currently a submodule within CMC, ensuring that any cluster simulations or binary populations are integrated with the same physics.


For Users
==================
.. toctree::
   :maxdepth: 2

   intro/index
   install/index
   examples/index
   output/index
   inifile/index


For Developers 
==================
.. toctree::
   :maxdepth: 2

   parallel/index
   src/index

