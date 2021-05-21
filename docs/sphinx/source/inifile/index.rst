.. _inifile:

###########################
Configuration files for CMC
###########################

The `cmc` command-line executable cannot run without a configuration file.  These ini files are designed to contain all the stellar evolution physics to generate a COSMIC population (i.e. the initial conditions for CMC) and all the dynamical options.  Note that there are two distinct modes you might want to run CMC in: either as a purely dynamics code, with no real stars or stellar evolution, or with all the bells and whistles, to produce realistic models of clusters.  These are sufficiently different that we have included two ini files in the examples directory:

    **KingProfile.ini** is designed to run a cluster with most of the available physics.  Use this for realistic star cluster simulations (it is what is described as "default" below).  The ini file can be found  `Here <https://github.com/ClusterMonteCarlo/CMC-COSMIC/blob/master/examples/KingProfile.ini>`_

    **PlummerSphere.ini** is designed to run a cluster with pure dynamics; no stars or stellar evolution, and will halt once the cluster reaches core collapse.  Use this for classical dynamics studies and learning the code.  The ini file can be found  `Here <https://github.com/ClusterMonteCarlo/CMC-COSMIC/blob/master/examples/PlummerSphere.ini>`_


.. warning::

    The parameter names are case insensitive, and we've tried to include everything currently used in the code here.  If you forget one, CMC will substitute a default value.  However, note that there are several default values still in the code that are NOT described here.  These are largely from either deprecated features, or parts of the code that are not working yet in the parallel version.  If they work at all, we would strongly advise against using them for any type of production run.  


[cmc]
=====


Although this is all one section, we have grouped the flags/parameters which get passed to the dynamics and binary stellar evolution codes into (rough) categories. Each group will start with a note to indicate. the type of parameter or flag.


.. _snapshotting:


SNAPSHOT FLAGS
--------------

===============================  =====================================================
``SNAPSHOTTING``                 Turn snapshotting on or of (doesn't control black hole snapshots below)

                                     ``0`` : Off

                                     ``1`` : On

                                 **SNAPSHOTTING = 1**

``SNAPSHOT_DELTACOUNT``          How many timesteps to write a snapshot

                                 **SNAPSHOT_DELTACOUNT = 1000**

``BH_SNAPSHOTTING``              Turn on snapshots for just the black holes (in the output.blackhole.snapshots.h5 file)

                                     ``0`` : Off

                                     ``1`` : On

                                 **BH_SNAPSHOTTING = 0**

``BH_SNAPSHOT_DELTACOUNT``       write the black hole snapshots every X timesteps

                                 **BH_SNAPSHOT_DELTACOUNT = 500**

``SNAPSHOT_CORE_COLLAPSE``       write additional snapshots around core collapse (in the output.snapshots.h5 file)

                                     ``0`` : Off

                                     ``1`` : On

                                 **SNAPSHOT_CORE_COLLAPSE = 0**

``SNAPSHOT_WINDOWS``             write snapshots every X units of time into the <output>.window.snapshots.h5 file.  Format is:

                                     start_w0,step_w0,end_w0;start_w1,step_w1,stop_w1 ... etc

                                 **SNAPSHOT_WINDOWS = "0,0.1,13.8"** (will write one snapshot to snapshot.h5 every 100 Myr)
 
``SNAPSHOT_WINDOW_UNITS``        What units should the snapshot windows be in.  Options are:

                                     ``Gyr`` : gigayears

                                     ``Trel`` : relaxation times

                                     ``Tcr`` : crossing times

                                 **SNAPSHOT_WINDOW_UNITS = Gyr**

===============================  =====================================================


RELAXATION FLAGS
----------------

===============================  =====================================================
``RELAXATION``                   perform two-body relaxation (0=off, 1=on)

                                     ``0`` : Off

                                     ``1`` : On

                                 **RELAXATION = 1**

``CMC_GAMMA``                    Value of the Coulomb Logarithm (lambda = log(gamma*N)).  Suggested values are:

                                     ``0.11`` : equal-mass clusters `Giersz & Heggie (1994) <https://ui.adsabs.harvard.edu/abs/1994MNRAS.268..257G/abstract>`_

                                     ``0.01`` : realistic initial mass functions `Rodriguez et al., (2018) <https://ui.adsabs.harvard.edu/abs/2018ComAC...5....5R/abstract>`_


                                 **CMC_GAMMA = 0.01**


``THETASEMAX``                   Maximum value of theta during two-body encounters.  Suggested values are:

                                     ``1.412`` : sqrt(2), for equal mass clusters

                                     ``1`` : for point-mass clusters


                                 **THETASEMAX = 1.412**
===============================  =====================================================


DYNAMICS FLAGS
--------------

==================================  =====================================================
``BINBIN``                          Turn on Fewbody encounters between two binaries 

                                     ``0`` : Off

                                     ``1`` : On

                                    **BINBIN = 1**


``BINSINGLE``                       Turn on Fewbody encounters between binaries and single objects

                                     ``0`` : Off

                                     ``1`` : On

                                    **BINSINGLE = 1**


``BH_CAPTURE``                      Turn on post-Newtonian corrections for black holes. NOTE: this activates GW captures for both single BHs and during fewbody encounters.  Note if SS_COLLISION=0 and BH_CAPTURE=1, GW captures will happen only in fewbody.

                                     ``0`` : Off

                                     ``1`` : On

                                    **BH_CAPTURE = 1**

``BH_RADIUS_MULTIPLYER``            Factor to multiply the radii of BHs by for collisions in fewbody (default is 5, since PN breaks down at ~10M)

                                    **BH_RADIUS_MULTIPLYER = 5**


``THREEBODYBINARIES``               Turn on three-body binary formation semi-analytic treatment from Morscher et al., 2013

                                     ``0`` : Off

                                     ``1`` : On

                                    **THREEBODYBINARIES = 1**



``ONLY_FORM_BH_THREEBODYBINARIES``  Only form three-body binaries if the three objects are black holes

                                     ``0`` : Off

                                     ``1`` : On

                                    **ONLY_FORM_BH_THREEBODYBINARIES = 1**


``MIN_BINARY_HARDNESS``             Minimum hardness ratio for forming three-body binaries (or breaking wide binaries if BINARY_BREAKING_MIN = 1)


                                    **MIN_BINARY_HARDNESS = 5**


``BINARY_BREAKING_MIN``             whether to use 10% of the interparticle seperation (0, default) or MIN_BINARY_HARDNESS (1) as the criterion for breaking wide binaries

                                     ``0`` : interparticle seperation

                                     ``1`` : MIN_BINARY_HARDNESS

                                    **BINARY_BREAKING_MIN = 0**


``SS_COLLISION``                    enable collisions between stars. NOTE: this activates collisions between single stars AND during fewbody encounters

                                     ``0`` : Off

                                     ``1`` : On

                                    **SS_COLLISION = 1**
 
``TIDAL_CAPTURE``                   Enable tidal captures bewteen individual stars.  Uses cross sections from Kim & Lee 1999 and Lombardi et al., 2006.  Only activated if SS_COLLISION = 1

                                     ``0`` : Off

                                     ``1`` : On

                                    **TIDAL_CAPTURE = 0**
 
``BHNS_TDE``                        Treat BH(NS)--MS TDEs in TDE vs direct collision limit.  Follows prescription in Kremer et al., 2020


                                     ``0`` : collision

                                     ``1`` : TDE

                                    **BHNS_TDE = 0**
==================================  =====================================================




INPUT OPTIONS
-------------

===============================  =====================================================
``INPUT_FILE``                   the input hdf5 file for our intial conditions generated using COSMIC

                                 **INPUT_FILE = input.hdf5**


===============================  =====================================================


TIDAL FIELD OPTIONS
-------------------

===============================  =====================================================
``TIDALLY_STRIP_STARS``          Tidally strip stars that pass the tidal boundary

                                     ``0`` : Off

                                     ``1`` : On

                                 **TIDALLY_STRIP_STARS = 1**

``TIDAL_TREATMENT``              choose the tidal cut-off criteria for removing stars

                                     ``0`` : radial criterion, :math:`{r_{\rm apo} > r_{t}}`

                                     ``1`` : Energy criterion, :math:`{\alpha E > \Phi(r_t)}`, following `Giersz et al. (2008) <https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..429G/abstract>`_

                                 **TIDAL_TREATMENT = 0**

``USE_TT_FILE``                  whether to take the tidal boundary from a file (i.e. a tidal tensor)

                                     ``0`` : Off

                                     ``1`` : On

                                 **USE_TT_FILE = 0**
 

``TT_FILE``                      name of tidal tensor file to take tidal boundary from.  Should be formatted as ``Time [Myr] T_xx T_yy T_zz T_xy T_xz T_yz [1/Myr^2]``

                                 **TT_FILE=NULL**
===============================  =====================================================


TERMINATION OPTIONS
-------------------

================================  =====================================================
``T_MAX_PHYS``                    maximum integration time (Gyr)

                                  **T_MAX_PHYS = 13.8**


``T_MAX``                         maximum integration time (relaxation time)

                                  **T_MAX = 100**


``T_MAX_COUNT``                   maximum number of timesteps

                                  **T_MAX_COUNT = 10000000**


``MAX_WCLOCK_TIME``               maximum amount of wallclock time to integrate for (seconds)

                                  **MAX_WCLOCK_TIME = 604800**


``CHECKPOINT_INTERVAL``           how often to save checkpoint files (seconds)

                                  **CHECKPOINT_INTERVAL = 7200**


``CHECKPOINTS_TO_KEEP``           how many previous checkpoints to keep at once

                                  **CHECKPOINTS_TO_KEEP = 2**


``TERMINAL_ENERGY_DISPLACEMENT``  energy change calculation stopping criterion (i.e. if :math:`e/e_0` changes by this much, stop calculation)

                                  **TERMINAL_ENERGY_DISPLACEMENT = 10**


``STOPATCORECOLLAPSE``            stop once the cluster reaches core collapse (< 100 stars in core)

                                     ``0`` : Off

                                     ``1`` : On

                                  **STOPATCORECOLLAPSE = 0**

``USE_DF_CUTOFF``                 whether to compute a lifetime according to dynamical friction criterion. This controls whether an external file (``DF_FILE``) is used

                                     ``0`` : Off

                                     ``1`` : On

                                  **USE_DF_CUTOFF = 0**

``DF_FILE``                       name of dynamical friction file to use as potential termination criterion.  Should be formatted as ``Time[Myr]`` ``Radius[Kpc]`` ``V_circ[km/s]`` ``M_enc[solMass]`` ``sigma[km/s]`` ``J[Kpc*km/s]``


                                  **DF_FILE = NULL**

``DF_INTEGRATED_CRITERION``       dynamical friction criterion to use for terminating the simulation

                                     ``0`` : :math:`{t_{\rm df} > \rm{TotalTime}}`

                                     ``1`` : :math:`{\int\frac{dt}{t_{\rm df}} > 1}`

                                  **DF_INTEGRATED_CRITERION = 1**

================================  =====================================================

OUTPUT OPTIONS
--------------

===============================  =====================================================
``MASS_PC``                      mass fractions for Lagrange radii, i.e. what fractions to actually print in the laggard files

                                 **MASS_PC = 0.0001,0.0003,0.0005,0.0007,0.0009,**
                                 **0.001,0.003,0.005,0.007,0.009,0.01,**
                                 **0.03,0.05,0.07,0.09,0.1,0.2,0.3,0.4,**
                                 **0.5,0.6,0.7,0.8,0.9,0.99**


``MASS_BINS``                    mass ranges for calculating derived quantities, i.e. what bins to use in mass for the different mass files

                                 **MASS_BINS = 0.1,1.0,10.0,100.0,1000.0**


``WRITE_EXTRA_CORE_INFO``        Write out information about cores that are defined differently from the standard

                                    ``0`` : Off

                                    ``1`` : On

                                 **WRITE_EXTRA_CORE_INFO = 0**


``WRITE_BH_INFO``                Write out information about BHs each timestep

                                    ``0`` : Off

                                    ``1`` : On

                                 **WRITE_BH_INFO = 1**


``WRITE_PULSAR_INFO``            Write out information about pulsars

                                    ``0`` : Off

                                    ``1`` : On

                                 **WRITE_PULSAR_INFO = 0**

``WRITE_MOREPULSAR_INFO``        Write a ton more information about neutron stars every PULSAR_DELTACOUNT timesteps


                                    ``0`` : Off

                                    ``1`` : On

                                 **WRITE_MOREPULSAR_INFO = 0**


``PULSAR_DELTACOUNT``            Pulsar output interval in time steps


                                 **PULSAR_DELTACOUNT = 1000**


``WRITE_STELLAR_INFO``           Write out information about stellar evolution for each single and binary star.  Warning, this creates a TON of output


                                    ``0`` : Off

                                    ``1`` : On

                                 **WRITE_STELLAR_INFO = 0**

===============================  =====================================================

CMC Parameters
--------------
Parameters for controlling the CMC run.  

.. warning::

    Don't touch these unless you know what you're doing.  And even then I probably wouldn't.  

===============================  =====================================================
``AVEKERNEL``                    one half the number of stars over which to average certain quantities

                                 **AVEKERNEL = 20**


``BH_AVEKERNEL``                 Same, but for three-body binary formation (which is fundamentally more local)

                                 **BH_AVEKERNEL = 3**


``MIN_CHUNK_SIZE``               minimum size of chunks that get partitioned across processors in the parallel code

                                 **MIN_CHUNK_SIZE = 40** (this is just twice the AVEKERNEL)

``IDUM``                         Random number generator seed.  Note this is different from the BSE RNG seed

                                 **IDUM = 1234**

``BSE_IDUM``                     The random number seed used by kick.f and setting initial pulsar spin period and magnetic field

                                 **BSE_IDUM = 1234**

``TIMER``                        Do profiling of the code, and print it out to the timers file.  Note that this introduces many MPI barriers

                                    ``0`` : Off

                                    ``1`` : On

                                 **TIMER = 0**              

``FORCE_RLX_STEP``               Force a relaxation timestep (useful when RELAXATION=0) 

                                    ``0`` : Off

                                    ``1`` : On

                                 **FORCE_RLX_STEP = 0**              

``DT_HARD_BINARIES``             calculate the binary interaction time steps by only considering hard binaries 

                                    ``0`` : Off

                                    ``1`` : On

                                 **DT_HARD_BINARIES = 0**              

``HARD_BINARY_KT ``              The minimum binary binding energy (in units of kT) for a binary to be considered 'hard' for the time step calculation.

                                 **HARD_BINARY_KT = 0.7**

``SAMPLESIZE``                   Number of samples keys to use for the parallel sample-sort algorithm

                                 **SAMPLESIZE = 1024**

``NUM_CENTRAL_STARS``            The number of central stars to use for calculating different qualities related to the timestep

                                 **NUM_CENTRAL_STARS = 300**
              

===============================  =====================================================


[bse]
=====

.. note::

    Although this is all one section, we have grouped the
    flags/parameters which get passed to the binary stellar evolution
    code into types. Each group will start with a note to indicate
    the type of parameter or flag.


SAMPLING FLAGS
--------------

=======================  =====================================================
``pts1``                 determines the timesteps chosen in each evolution phase as
                         decimal fractions of the time taken in that phase for
                         Main Sequence (MS) stars

                         **pts1 = 0.001** following `Bannerjee+2019 <https://ui.adsabs.harvard.edu/abs/2019arXiv190207718B/abstract>`_

``pts2``                 determines the timesteps chosen in each evolution phase as
                         decimal fractions of the time taken in that phase for
                         Giant Branch (GB, CHeB, AGB, HeGB) stars

                         **pts2 = 0.01** following `Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_

``pts3``                 determines the timesteps chosen in each evolution phase as
                         decimal fractions of the time taken in that phase for
                         HG, HeMS stars

                         **pts3 = 0.02** following `Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_
=======================  =====================================================



METALLICITY FLAGS
-----------------

=======================  =====================================================
``zsun``                 Sets the metallicity of the Sun which primarily affects
                         stellar winds.

                         **zsun = 0.014** following `Asplund 2009 <https://ui.adsabs.harvard.edu/abs/2009ARA%26A..47..481A/abstract>`_
=======================  =====================================================


WIND FLAGS
----------

=======================  =====================================================
``windflag``             Selects the model for wind mass loss for each star

                            ``0`` : Standard SSE/BSE (`Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_)

                            ``1`` : StarTrack (`Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_)

                            ``2`` : Metallicity dependence for O/B stars and Wolf Rayet stars (`Vink+2001 <http://adsabs.harvard.edu/abs/2001A&amp;A...369..574V>`_, `Vink+2005 <https://ui.adsabs.harvard.edu/abs/2005A%26A...442..587V/abstract>`_)

                            ``3`` : Same as 2, but LBV-like mass loss for giants
                            and non-degenerate stars beyond the
                            Humphreys-Davidson limit

                         **windflag = 3**

``eddlimflag``           Limits the mass-loss rate of low-metallicity stars near
                         the Eddington limit
                         (see `Grafener+2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...535A..56G/abstract>`_, `Giacobbo+2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.2959G/abstract>`_).

                            ``0`` : does not apply Eddington limit

                            ``1`` : applies Eddington limit

                         **eddlimflag = 0**

``neta``                 Reimers mass-loss coefficent (`Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_).
                         Note: this equation has a typo. There is an extra
                         :math:`{\eta}` out front; the correct rate is directly proportional
                         to :math:`{\eta}`.
                         See also `Kurdritzki+1978, Section Vb <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1978A%26A....70..227K&link_type=ARTICLE&db_key=AST&high=#page=12>`_ for discussion.

                            ``positive value`` : supplies :math:`{\eta}` to `Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_

                         **neta = 0.5**

``bwind``                Binary enhanced mass loss parameter.
                         See `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                            ``positive value`` : supplies B\ :sub:`w` to `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_

                         **bwind = 0, inactive for single**

``hewind``               Helium star mass loss parameter: 10\ :sup:`-13` *hewind* L\ :sup:`2/3` gives He star mass-loss. Equivalent to 1 - :math:`{\mu}` in the last equation on `page 19 of SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_.

                         **hewind = 0.5**

``beta``                 Wind velocity factor: v\ :sub:`wind` :sup:`2` goes like *beta*. See `Equation 9 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                            ``negative value`` : StarTrack (`Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_)

                            ``positive value`` : supplies :math:`{\beta}`\ :sub:`w` to `Equation 9 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_

                         **beta = -1.0**

``xi``                   Wind accretion efficiency factor, which gives the fraction
                         of angular momentum lost via winds from the primary that
                         transfers to the spin angular momentum of the companion.
                         Corresponds to :math:`{\mu}`\ :sub:`w` in `Equation 11 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                            ``positive value`` : supplies :math:`{\mu}`\ :sub:`w` in `Equation 11 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_

                         **xi = 0.5**

``acc2``                 Bondi-Hoyle wind accretion factor where the mean wind accretion rate onto the secondary is proportional to *acc2*. See `Equation 6 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_.

                            ``positive value`` : supplies :math:`{\alpha}`\ :sub:`w` in `Equation 6 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_

                         **acc2 = 1.5**
=======================  =====================================================


COMMON ENVELOPE FLAGS
---------------------

**Note:** there are cases where a common envelope is forced regardless of the
critical mass ratio for unstable mass transfer. In the following cases, a
common envelope occurs regardless of the choices below:

**contact** : the stellar radii go into contact (common for similar ZAMS systems)

**periapse contact** : the periapse distance is smaller than either of the stellar radii (common for highly eccentric systems)

**core Roche overflow** : either of the stellar radii overflow their component's Roche radius (in this case, mass transfer from the convective core is always dynamically unstable)

=======================  =====================================================
``alpha1``               Common-envelope efficiency parameter which scales the
                         efficiency of transferring orbital energy to the
                         envelope. See `Equation 71 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

                            ``positive values`` : supplies :math:`{\alpha}` to `Equation 71 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_

                         **alpha1 = 1.0**

``lambdaf``              Binding energy factor for common envelope evolution.
                         The initial binding energy of the stellar envelope
                         goes like 1 / :math:`{\lambda}`. See `Equation 69 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

                            ``positive values`` : uses variable lambda prescription detailed
                            in appendix of `Claeys+2014 <https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract>`_
                            where lambdaf is the fraction of the ionization energy that can go into ejecting
                            the envelope; to use this prescription without extra ionization energy, set lambdaf=0

                            ``negative values`` : fixes :math:`{\lambda}` to a value of -1.0* *lambdaf*

                         **lambdaf = 0.0**

``ceflag``               Selects the `de Kool 1990 <https://ui.adsabs.harvard.edu/abs/1990ApJ...358..189D/abstract>`_
                         model to set the initial orbital energy using the
                         total mass of the stars instead of the core masses as
                         in `Equation 70 of Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_.

                            ``0`` : Uses the core mass to calculate initial
                            orbital energy as
                            in `Equation 70 of Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``1`` : Uses the `de Kool 1990 <https://ui.adsabs.harvard.edu/abs/1990ApJ...358..189D/abstract>`_
                            model

                         **ceflag = 0**

``cekickflag``           Selects which mass and separation values to use when
                         a supernova occurs during the CE and a kick
                         needs to be applied.

                            ``0`` : uses pre-CE mass and post-CE sep (BSE default)

                            ``1`` : uses pre-CE mass and sep values

                            ``2`` : uses post-CE mass and sep

                         **cekickflag = 2**

``cemergeflag``          Determines whether stars that begin a CE
                         without a core-envelope boundary automatically lead to
                         merger in CE. These systems include:
                         kstars = [0,1,2,7,8,10,11,12].

                            ``0`` : allows the CE to proceed

                            ``1`` : causes these systems to merge in the CE

                         **cemergeflag = 0**

``cehestarflag``         Uses fitting formulae from `Tauris+2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.2123T/abstract>`_
                         for evolving RLO systems with a helium star donor
                         and compact object accretor.
                         NOTE: this flag will override choice made by
                         cekickflag if set

                            ``0`` : does NOT use Tauris+2015 at all

                            ``1`` : uses Tauris+2015 fits for final period only

                            ``2`` : uses Tauris+2015 fits for both final mass and final period

                         **cehestarflag = 0**

``qcflag``               Selects model to determine critical mass ratios for the
                         onset of unstable mass transfer and/or a common envelope
                         during RLO.
                         NOTE: this is overridden by qcrit_array if any of the
                         values are non-zero.

                            ``0`` : follows `Section 2.6 of Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                            (Default BSE)

                            ``1`` : same as 0 but with `Hjellming & Webbink 1987 <https://ui.adsabs.harvard.edu/abs/1987ApJ...318..794H/abstract>`_
                            for GB/AGB stars

                            ``2`` : follows `Table 2 of Claeys+2014 <https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract>`_

                            ``3`` : same as 2 but with `Hjellming & Webbink 1987 <https://ui.adsabs.harvard.edu/abs/1987ApJ...318..794H/abstract>`_
                            for GB/AGB stars

                            ``4`` : follows `Section 5.1 of Belcyznski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_ except for WD donors which follow BSE

                            ``5`` : follows `Section 2.3 of Neijssel+2020 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract>`_ Mass transfer from stripped stars is always assumed to be dynamically stable

                         **qcflag = 1**

                         .. csv-table:: Comparison of Q Crit Values (Donor Mass/Accretor Mass) For Each Donor Kstar Type Across Flag Options
                            :file: qcrit_table.csv
                            :header-rows: 1


                         Eq.1: ``qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))``, which is from Hjellming & Webbink 1983

                         Eq.2: ``qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0``, which is from Claeys+ 2014

``qcrit_array``          Array with length: 16 for user-input values for the
                         critical mass ratios that govern the onset of unstable
                         mass transfer and a common envelope. Each item is set
                         individually for its associated kstar, and a value of
                         0.0 will apply prescription of the qcflag for that kstar.

                         **Note:** there are cases where a common envelope is forced
                         regardless of the critical mass ratio for unstable mass
                         transfer; in the following cases, a common envelope occurs
                         regardless of the qcrit or qcflag

                         **qcrit_array = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]**

=======================  =====================================================


KICK FLAGS
----------

=======================  =====================================================
``kickflag``             Sets the particular natal kick prescription to use
                         Note that ``sigmadiv``, ``bhflag``, ``bhsigmafrac``,
                         ``aic``, and ``ussn``, which are described below, are
                         only used when ``kickflag=0``

                            ``0`` : The standard COSMIC kick prescription, where
                            kicks are drawn from a bimodal distribution with
                            standard FeCCSN getting a kick drawn from a Maxwellian
                            distribution with dispersion parameter ``sigma`` and ECSN
                            are drawn according to ``sigmadiv``. This setting has
                            additional possible options for ``bhflag``, ``bhsigmafrac``,
                            ``aic`` and ``ussn``.

                            ``-1`` : Natal kicks are drawn according to ``sigma`` and
                            scaled by the ejecta mass and remnant mass following Eq. 1 of
                            `Giacobbo & Mapelli 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract>`_

                            ``-2`` : Natal kicks are drawn according to ``sigma`` and
                            scaled by just the ejecta mass following Eq. 2 of
                            `Giacobbo & Mapelli 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract>`_

                            ``-3`` : Natal kicks are drawn according to Eq. 1 of
                            `Bray & Eldridge 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.461.3747B/abstract>`_

                         **default=0**

``sigma``                Sets the dispersion in the Maxwellian for the
                         SN kick velocity in km/s

                            ``positive value`` : sets Maxwellian dispersion

                         **default=265.0**

``bhflag``               Sets the model for how SN kicks are applied to BHs
                         where bhflag != 0 allows velocity kick at BH formation

                            ``0`` : no BH kicks

                            ``1`` : fallback-modulated kicks following
                            `Fryer+2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`_

                            ``2`` : kicks decreased by ratio of BH mass to NS mass
                            (1.44 Msun); conserves linear momentum

                            ``3`` : full strength kick drawn from Maxwellian
                            with dispersion = *sigma* selected above

                         **bhflag = 1**

``ecsn``                 Allows for electron capture SN and sets the
                         maximum ECSN mass range at the time of SN

                            ``0`` : turns off ECSN

                            ``positive values`` : `BSE (Hurley+2002) <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                            and `StarTrack (Belczynski+2008) <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_
                            use ecsn = 2.25, while `Podsiadlowksi+2004 <https://ui.adsabs.harvard.edu/abs/2004ApJ...612.1044P/abstract>`_
                            use ecsn = 2.5

                         **ecsn = 2.5**

``ecsn_mlow``            Sets the low end of the ECSN mass range

                            ``positive values`` : `BSE (Hurley+2002) <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                            use ecsn_mlow = 1.6, while `StarTrack (Belczynski+2008) <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_
                            use ecsn_mlow = 1.85, while `Podsiadlowksi+2004 <https://ui.adsabs.harvard.edu/abs/2004ApJ...612.1044P/abstract>`_
                            use ecsn_mlow = 1.4

                         **ecsn_mlow = 1.4**

``sigmadiv``             Sets the modified ECSN kick strength

                         ``positive values`` : divide *sigma* above by *sigmadiv*

                         ``negative values`` : sets the ECSN *sigma* value

                         **sigmadiv = -20.0**

``aic``                  reduces kick strengths for accretion induced collapse SN
                         according to *sigmadiv*

                            ``0`` : AIC SN receive kicks drawn from Maxwellian
                            with dispersion = *sigma* above

                            ``1`` : sets kick strength according to *sigmadiv*
                            NOTE: this will applies even if ecsn = 0.0

                         **aic = 1**

``ussn``                 Reduces kicks according to the *sigmadiv* selection
                         for ultra-stripped supernovae which happen whenever
                         a He-star undergoes a CE with a compact companion

                            ``0`` : USSN receive kicks drawn from Maxwellian
                            with dispersion = *sigma* above

                            ``1`` : sets kick strength according to *sigmadiv*

                         **ussn = 0**

``pisn``                 Allows for (pulsational) pair instability supernovae
                         and sets either the model to use or the maximum mass
                         of the remnant.

                            ``0`` : no pulsational pair instability SN

                            ``-1`` : uses the formulae from `Spera & Mapelli 2017 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.4739S/abstract>`_

                            ``-2`` : uses a polynomial fit to `Table 1 in Marchant+2018 <https://ui.adsabs.harvard.edu/abs/2018arXiv181013412M/abstract>`_

                            ``-3`` : uses a polynomial fit to `Table 5 in Woosley 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...878...49W/abstract>`_

                            ``positive values`` : turns on pulsational pair
                            instability SN and sets the maximum mass of the allowed
                            remnant

                         **pisn = 45.0**

``bhsigmafrac``          Sets a fractional modification which scales down *sigma*
                         for BHs. This works in addition to whatever is chosen for
                         *bhflag*, and is applied to *sigma* **before** the *bhflag*
                         prescriptions are applied

                            ``values between [0, 1]`` : reduces *sigma* by *bhsigmafrac*

                         **bhsigmafrac = 1.0**

``polar_kick_angle``     Sets the opening angle of the SN kick relative to the
                         pole of the exploding star, where 0 gives strictly polar
                         kicks and 90 gives fully isotropic kicks

                            ``values between [0, 90]`` : sets opening angle for SN kick

                         **polar_kick_angle = 90.0**

``natal_kick_array``     Array of dimensions: (2,5) which takes user input values
                         for the SN natal kick, where the first row corresponds to the
                         first star and the second row corresponds to the second star and
                         columns are: [vk, phi, theta, mean_anomaly, rand_seed].
                         NOTE: any numbers outside these ranges will be sampled
                         in the standard ways detailed above.

                            ``vk`` : valid on the range [0, inf]

                            ``phi`` : co-lateral polar angle in degrees, valid from
                            [-90, 90]

                            ``theta`` : azimuthal angle in degrees, valid from
                            [0, 360]

                            ``mean_anomaly`` : mean anomaly in degrees,
                            valid from [0, 360]

                            ``rand_seed`` : supplied if restarting evolution after
                            a supernova has already occurred

                         **natal_kick_array = [[-100.0,-100.0,-100.0,-100.0,0.0][-100.0,-100.0,-100.0,-100.0,0.0]]**
=======================  =====================================================

    
REMNANT MASS FLAGS
------------------

===================  =====================================================
``remnantflag``      Determines the remnant mass prescription used for NSs and BHs.

                            ``0`` : follows `Section 6 of Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_
                            (default BSE)

                            ``1`` : follows `Belczynski+2002 <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..407B/abstract>`_

                            ``2`` : follows `Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                            ``3`` : follows the rapid prescription from `Fryer+2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`_, with updated proto-core mass from `Giacobbo & Mapelli 2020 <https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract>`_

                            ``4`` : delayed prescription from `Fryer+2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`_

                     **remnantflag = 3**

``mxns``             Sets the boundary between the maximum NS mass
                     and the minimum BH mass

                            ``positive values`` : sets the NS/BH mass bounary

                     **mxns = 3.0**

``rembar_massloss``  Determines the prescriptions for mass conversion from
                     baryonic to gravitational mass during the collapse of
                     the proto-compact object

                            ``positive values`` : sets the maximum amount of mass loss, which should be about 10% of the maximum mass of an iron core (:math:`{\sim 5 \mathrm{M}_\odot}` Fryer, private communication)

                            ``-1 < *rembar_massloss* < 0`` : assumes that proto-compact objects lose a constant fraction of their baryonic mass when collapsing to a black hole (e.g., *rembar_massloss* = -0.1 gives the black hole a gravitational mass that is 90% of the proto-compact object's baryonic mass)

                     **rembar_massloss = 0.5**
===================  =====================================================


REMNANT SPIN FLAGS
------------------

=======================  ===============================================================
``bhspinflag``           Uses different prescriptions for BH spin after formation

                            ``0`` : sets all BH spins to *bhspinmag*

                            ``1`` : draws a random BH spin between 0 and bhspinmag for every BH

                            ``2`` : core-mass dependent BH spin (based on `Belczynski+2017 v1 <https://arxiv.org/abs/1706.07053v1>`_)

                         **bhspinflag = 0**

``bhspinmag``            Sets either the spin of all BHs or the upper limit of the uniform distribution for BH spins

                            ``values >= 0.0`` : spin or upper limit value

                         **bhspinmag = 0.0**
=======================  ===============================================================


GR ORBITAL DECAY FLAG
---------------------

.. note::

    In CMC, GR orbital decay is handled separately from BSE for binary black holes, and is unaffected by the below flag


=======================  ===============================================================
``grflag``               Turns on or off orbital decay due to gravitational wave radiation

                            ``0`` : No orbital decay due to GR

                            ``1`` : Orbital decay due to GR is included

                         **grflag = 1**
=======================  ===============================================================


MASS TRANSFER FLAGS
-------------------

=======================  =====================================================
``eddfac``               Eddington limit factor for mass transfer.

                            ``1`` : mass transfer rate is limited by the
                            Eddington rate following Equation 67 in
                            `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``values >1`` : permit super-Eddington accretion
                            up to value of *eddfac*

                         **eddfac = 1.0**

``gamma``                Angular momentum prescriptions for mass lost during RLO
                         at super-Eddington mass transfer rates

                            ``-1`` : assumes the lost material carries away the
                            specific angular momentum of the primary

                            ``-2`` : assumes material is lost from the system as
                            if it is a wind from the secondary

                            ``>0`` : assumes that the lost material takes away a
                            fraction *gamma* of the orbital angular momentum

                         **gamma = -2.0**

``don_lim``              Calculates the rate of thermal mass loss through Roche
                         overflow mass transfer from the donor star

                            ``-1`` : donor mass loss rate is calculated following
                            `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``-2`` : donor mass loss rate is calculated following
                             `Claeys+2014 <https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract>`_

``acc_lim``              Limits the amount of mass accreted during Roche overflow

                            ``-1`` : limited to 10x's the thermal rate of the accretor
                            for MS/HG/CHeB and unlimited for GB/EAGB/AGB stars

                            ``-2`` : limited to 1x's the thermal rate of the accretor
                            for MS/HG/CHeB and unlimited for GB/EAGB/AGB stars

                            ``-3`` : limited to 10x's the thermal rate of the accretor
                            for all stars

                            ``-4`` : limited to 1x's the thermal rate of the accretor
                            for all stars

                            ``>=0`` : sets overall accretion fraction of donor mass
                            as in Belcyznski+2008 w/ acc_lim = 0.5

=======================  =====================================================


TIDES FLAGS
-----------

=======================  =====================================================
``tflag``                Activates tidal circularisation following
                         `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``0`` : no tidal circularization

                            ``1`` : activates tidal circularization

                         **tflag = 1**

``ST_tide``              Activates StarTrack setup for tides following
                         `Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                            ``0`` : follows `BSE <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``1`` : follows `StarTrack <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                         **ST_tide = 1**

``fprimc_array``         controls the scaling factor for convective tides
                         each item is set individually for its associated kstar
                         The releveant equation is `Equation 21 of Hurley+2002 <https://watermark.silverchair.com/329-4-897.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAnAwggJsBgkqhkiG9w0BBwagggJdMIICWQIBADCCAlIGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMYUoYtydpxVKmZePqAgEQgIICI1b5IZldHg9_rX6JacIe-IR042LnNi-4F9DMp-2lm3djjQ8xehKOv5I0VBjSNJfa6n-FErAH7ed1llADY7tMDTvqo1GHKBMDslNku5XDGfmae0sF-Zp5ndeGoZsyqISABLHEbdY4VFl8Uz_6jzAuBjGztnuxVmUh9bKIOaxuDpfB3Mn2xOfP9lcCVkjzQ0JWzr98nQNmVwDkI9bPv98Ab46BjBdGdcBKajCC-sqASjtmAQS2h6SGTTBqyRAyigqXcPtWf3Ye1SbxtL3zag6_Lf01rgCoUCK9eT_pavb5F8vVkUTMWbZQ79DWxn5pfZYi72C7_BtlPoUnS8Gs3wvw18BTIaHTKblwh225DcXuTEh_ngMmRvPEVctvG8tjlr9md-eFK0cEsq0734eGYtnwxeqvFxcWsW6mRbXrFHFsInQK16j6n36XuCimY665l_-HPAuu-lTTlwpMTUR7K1eYMBsco_tp_TdxEipRNvBpaWZX3J0FxPMzi84Y01UvWiW69pxb-LLTpf8aG4YCm9asRFyfDZ9nbSdgrIlCiuzy7QSmkvsHOaTEecmwRimFRycDuIuWLvA_tILmYCIM2KzvqYJSVCQPJH39xEHZG8LbMqImwAVYO3H90qh-90gNrtZn4ofSskcgqxeqfZly9CPfmEevX5s-SlLHMh1N6gdZwenvMC0kTWg_rskbvGiANtuGngD-kKDbunGpYJU_nI7uDnhGtdY#page=5>`_

                            ``positive values`` : sets scaling factor of
                            Equation 21 referenced above

                         **fprimc_array = [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                         2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                         2.0/21.0,2.0/21.0]**
=======================  =====================================================


WHITE DWARF FLAGS
-----------------

=======================  =====================================================
``ifflag``               Activates the initial-final white dwarf mass relation
                         from Han+1995 `Equations 3, 4, and 5 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1995MNRAS.272..800H&link_type=ARTICLE&db_key=AST&high=#page=4>`_.

                            ``0`` : no modifications to BSE

                            ``1`` : activates initial-final WD mass relation

                         **ifflag = 0**

``wdflag``               Activates an alternate cooling law found in the description
                         immediately following `Equation 1 <http://iopscience.iop.org/article/10.1086/374637/pdf#page=3>`_
                         in Hurley & Shara 2003.
                         Equation 1 gives the BSE default Mestel cooling law.

                            ``0`` : no modifications to BSE

                            ``1`` : activates modified cooling law

                         **wdflag = 1**

``epsnov``               Fraction of accreted matter retained in a nova eruption.
                         This is relevant for accretion onto degenerate objects;
                         see Section 2.6.6.2 in `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_.

                            ``positive values between [0, 1]`` : retains *epsnov*
                            fraction of accreted matter

                         **epsnov = 0.001**
=======================  =====================================================

PULSAR FLAGS
------------

=======================  =====================================================
``bdecayfac``            Activates different models for accretion induced field decay; see
                         `Kiel+2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

                            ``0`` : uses an exponential decay

                            ``1`` : uses an inverse decay

                         **bdecayfac = 1**

``bconst``               Sets the magnetic field decay time-scale for pulsars following
                         Section 3 of `Kiel+2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

                            ``negative values`` : sets k in Myr from Equation 8 to
                            -1 * *bconst*

                         **bconst = -3000**

``ck``                   Sets the magnetic field decay time-scale for pulsars following
                         Section 3 of `Kiel+2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

                            ``negative values`` : sets :math:`{\tau}`\ :sub:`b` in Myr
                            from Equation 2 to  -1 * *ck*

                         **ck = -1000**
=======================  =====================================================

MIXING VARIABLES
----------------

=======================  =====================================================

``rejuv_fac``            Sets the mixing factor in main sequence star collisions.
                         This is hard coded to 0.1 in the original BSE release
                         and in Equation 80 of `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                         but can lead to extended main sequence lifetimes in some cases.

                             ``positive values`` : sets the mixing factor

                         **rejuv_fac = 1.0**

``rejuvflag``            Sets whether to use the orginal prescription for mixing
                         of main-sequence stars (based on equation 80 of `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_)
                         or whether to use the ratio of the pre-merger He core
                         mass at the base of the giant branch to the merger product's
                         He core mass at the base of the giant branch


                            ``0`` : no modifications to BSE

                            ``1`` : modified mixing times

                         **rejuvflag = 0**

``bhms_coll_flag``       If set to 1 then if BH+star collision and if Mstar > Mbh, do not destroy the star

                         **default = 0**

=======================  =====================================================


MAGNETIC BRAKING FLAGS
----------------------

=======================  =====================================================
``htpmb``                Activates different models for magnetic braking

                            ``0`` : no modifications to BSE

                            ``1`` : follows `Ivanona and Taam 2003 <https://ui.adsabs.harvard.edu/abs/2003ApJ...599..516I/abstract>`_

                         **htpmb = 1**
=======================  =====================================================



MISC FLAGS
----------

=======================  =====================================================
``ST_cr``                Activates different convective vs radiative boundaries

                            ``0`` : no modifications to BSE

                            ``1`` : follows `StarTrack <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                         **ST_cr = 1**
=======================  =====================================================
