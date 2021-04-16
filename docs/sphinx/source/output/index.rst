.. _output:

####################
Analyzing CMC Output
####################



.. ipython::
   :suppress:

   In [143]: import os

   In [144]: from matplotlib.pylab import *

   In [145]: style.use(os.path.abspath('./source/matplotlib_style_file'))

   In [146]: ion() 


Here we list all of the output files that are generated in a typical CMC run, and offer some suggestsions on how to analyze them.  Generally, the output comes in one of two types: `log files`, where individual events (e.g. stellar collisions, mergers, fewbody interactions) are logged whenenver they happen or when specific quantities are recorded every timestep (e.g. the number of black holes, the core radius, the lagrange radii), and `snapshots`, where the state of every star and binary in the cluster is recorded.  These are saved in plain text dat files and hdf5 files, respectively.

Most of the log files are in plain text (.dat) with easily readable columns that can be loaded in python with either ``loadtxt`` or the appropriate pandas commands.  However, a handful of the event log files, such as the binary interaction (binint.log), collision (collision.log), and stellar evolution distruption (semergedistrupt.log), are somewhat more difficult to parse.  We have provided a python package to make parsing these easier; it is located in the ``CMC-COSMIC/tools`` directory of the main CMC folder `(the cmc_parser.py file) <https://github.com/ClusterMonteCarlo/CMC-COSMIC/tree/master/tools>`_.

The HDF5 snapshots are designed to be importable as Pandas tables and are described below.  These snapshots can also be converted into various observational quantities, such as surface brightness profiles, velocity dispersion profiles, and more, using the `cmctoolkit` package (`Rui et al., 2021 <https://ui.adsabs.harvard.edu/abs/2021arXiv210305033R/abstract>`_).  We include the python code from the main package in the ``CMC-COSMIC/tools`` folder, but please see `here <https://github.com/NicholasRui/cmctoolkit>`_ for complete instructions and useage.

.. note::

        Here we are using the output from the :ref:`realistic-initial-conditions` example, running it for 15 Myr (though note that the default ini file runs it for 13.8 Gyr) and having it print a window snapshot every 10 Myr.  On 28 cores of the Vera machine, this takes about 10 minutes.  

==========
Code Units
==========

CMC uses the standard N-body units `(Henon 1971a) <https://link.springer.com/article/10.1007/BF00649159>`_ in all numerical calculations: the gravitational constant :math:`{G=1}`, the initial total mass :math:`{M=1}`, and the initial total energy :math:`{E_0=-1/4}`. The conversion from code units to physical units for time is done by calculating the initial half-mass relaxation time in years. The file ``<output>.conv.sh`` contains all conversion factors used to convert from code to physical units. For example, to convert from code units of time to Myr, multiply by the factor timeunitsmyr, and so on. 

This file can be parsed using the ``cmc_parser.py`` file for easier use in python:

.. ipython:: python

        import cmc_parser as cp
        conv_file = cp.conversion_file('source/example_output/output.conv.sh')
        print(conv_file.mass_msun)
        print(conv_file.__dict__)

==========
Log Files
==========



initial.dyn.dat
---------------
The `initial.dyn.dat` files contains various theoretical quantities pertaining to the dynamical evolution of your cluster, such as theoretical core radius, total mass, central density, central velocity dispersion, etc. This file may be used to create, for example, plots of core radii or half-light radii versus time for a given simulation. 

================  =====================================================
``t``               time
``dt``              timestep
``tcount``
``N``
``M``
``VR``               Measure how far away the cluster is from virial equilibrium (-2.0*Etotal.K/Etotal.P)
``N_c``              core
``r_c``
``r_max``            max radius of star in cluster
``Etot``
``KE``
``PE``
``Etot_int``        Total internal energy of stars
``Etot_bin``
``E_cenma``         Central BH mass (i.e., when there is a central IMBH)
``Eesc``      escape energy
``Ebesc``
``Eintesc``         Total internal energy in the escaped stars
``Eoops``           Totoal excess energy from breaking up triples and soft binaries
``Etot+Eoops``
``r_h``
``rho_0``
``rc_spitzer``      Core radius as defined in Spitzer 1987(sqrt(3.0 * sqr(central.v_rms) / (4.0 * PI * central.rho)))
``v0_rms``          rms velo dispersion at center?
``rc_nb``           Core radius calculated with density weighted averages as in Casertano & Hut (1985)
``DMse``			     Total cmass loss from the cluster per time step[:math:`{M_{\odot}}`]
``DMrejuv`` 	     Mass loss from rejuvenation per time step[:math:`{M_{\odot}}`]
``N_c_nb``          Number of stars in the core (4.0 / 3.0 * PI * cub(rc_nb) * (central.n / 2.0))
================  =====================================================

initial.binint.log
------------------

Over the course of the evolution of the cluster, single stars and binaries will frequently undergo three- and four-body dynamical encounters, which are integrated directly in CMC using the Fewbody package (Fregeau et al. 2007). The file `initial.binint.log` records all input and output parameters (e.g., component masses, IDS, stellar types, semi-major axes, etc.) each type fewbody is called. 

Every encounter information is printed between two lines of asterisks.
Below is an exemplary output::

>>> ********************************************************************************
>>> type=BS t=5.85010072e-06
>>> params: b=1.46611 v=0.379587
>>> input: type=single m=0.0284732 R=0.215538 Eint=0 id=170307 ktype=0
>>> input: type=binary m0=0.211113 m1=0.148022 R0=0.22897 R1=0.170412 Eint1=0 Eint2=0 
>>> ... id0=33128 id1=1255329 a=0.0908923 e=0.0641548 ktype1=0 ktype2=0
>>> status: DE/E=-1.79889e-08 DE=1.71461e-10 DL/L=2.54957e-08 DL=8.18406e-10 DE_GW/E=-0
>>> ... DE_GW=0 v_esc_cluster[km/s]=77.9847 tcpu=0.01
>>> outcome: nstar=3 nobj=2:  0 [1 2] (single-binary)
>>> output: type=single m=0.0284732 R=0.215538 Eint=0 id=170307 ktype=0
>>> output: type=binary m0=0.211113 m1=0.148022 R0=0.22897 R1=0.170412 Eint1=0 Eint2=0
>>> ... id0=33128 id1=1255329 a=0.09094 e=0.123848 ktype1=0 ktype2=0
>>> ********************************************************************************

==============================  =====================================================
``type``						encounter type (BS or BB)
``t``							encounter time
``b``							impact parameter
``v``							relative velocity
``m``							mass [:math:`{M_{\odot}}`]
``R``							radius
``Eint``			            
``id``							ID number
``kytpe``						stellar type
``a``							semi-major axis
``e``							eccentricity
``dE/E``			
``DE``
``DL/L``
``DL``
``DE_GW``
``v_esc_cluster``				escape velocity [km/s]
``tcpu``
``nstar``						number of stars
``nobj``						number of objects (single/binary)
``i [j k]``					    final configuration after encounter, e.g.,  0 [1 2] (single-binary)
==============================  =====================================================

Objects are labelled starting from 0 to 3. The binary-single and binary-binary encounters are denoted as BS and BB, respectively. For type=binary, indices 0 and 1 in mass, radius,id,etc. denote the primary and secondary objects in a binary.

Possible outcomes for ``type=BS``:

* single-binary 0 [1 2]
* binary-single [2 0] 1
* single-single-single 0 1 2
* single-single 0:1 2
* binary [0:1 2]
* single 0:1:2

Possible outcomes for ``type=BB``: 

* binary [0 1:2:3]
* single-binary 0:1 [2 3]
* binary-single [0:1 3] 2
* binary-binary [0 1] [2 3] 
* single-triple 0 [[1 3] 2]
* triple-single [[0 1] 3] 2
* single-single-binary 3 1 [2 0]
* binary-single-single [0 1] 3 2
* single-binary-single 0 [1 3] 2

0:1 denotes the fact that objects 0 and 1 have merged, and [0 1] indicates that objects 0 and 1 have formed a binary. The same is true for any pairs from 0 to 3.


initial.bh.dat
--------------

This file contains the number of BHs (as well as BH binaries, etc.) at each dynamical time step. This is useful to plot, for example, the number of retained BHs versus time. 

==============================  =====================================================
``tcount``						 time count
``Totaltime``					 total time
``Nbh,tot``						 total # BH
``Nbh,single``					 # single BH
``Nbinarybh``					 # binary BH
``Nbh-bh``						 # BH-BH binaries
``Nbh-nonbh``			         # BH-non BH binaries
``Nbh-ns``					     # BH-NS binaries
``Nbh-wd``						 # BH-WD binaries
``N_bh-star``				 	  including MS, giants not compact objs
``Nbh-ms``						 # BH-MS binaries	
``Nbh-postms``			         # BH-postMS binaries
``fb_bh``						 # binaries containing a black hole / total # systems containing a black hole
==============================  =====================================================



initial.BBHmerger.dat
---------------------

#List of all BBH mergers
  2 #0)model_num 1)rv 2)rg 3)Z 4)N 5)merger_time(Myr) 6)id1 7)id2 8)m1 9)m2 1    0)merger_channel
  3 #Merger channels: 1)Ejected 2)In-cluster 2-body 3)In-cluster 3-body 4)In-    cluster 4-body 5)In-cluster single-single capture

This file lists masses and merger times for all binary BH mergers occurring in the given simulation. Also listed are merger channels (ejected, in-cluster, etc.) as described in Section 9 of Kremer et al. 2019. 


initial.collision.log
---------------------

This file lists stellar types and properties for all stellar collisions occurring in the given simulation. See Sections 6 and 7 of Kremer et al. 2019 for further detail. 

==============================  =====================================================
``t``						     collision time
``interaction type``		     e.g., single-binary, binary-binary, etc.
``idm(mm)``						 ID_merger(mass of the merged body)
``id1(m1)``					     ID_1 (mass of the collided body_1)
``id2(m2)``					 	 ID_2 (mass of the collided body_2)
``r``						 
``typem``			             merger stellar type
``type1``					     body_1 stellar type
``type2``						 body_2 stellar type
``b``                            [:math:`R_{\odot}`]
``vinf``                         [km/s]
``rad1``                         radius of body_1
``rad2``                         radius of body_2
``rperi``                        pericenter distance at collision
``coll_mult``                     colliison multipole if sticky sphere = 1 if bh ns TDE (>1)
==============================  =====================================================


The single-single, binary-single, etc tells whether the collision occurred during a binary encounter or not. When there are three stars listed for the collision it means that all three stars collided during the encounter. This is rare, but it does happen occasionally. Typically, one will see something like::

>>> t=0.00266079 binary-single idm=717258(mm=1.0954) id1=286760(m1=0.669391):id2=415309 
>>> (m2=0.426012) (r=0.370419) typem=1 type1=0 type2=0

In this case the colliding stars are m1=0.66 and m2=0.42. The information about the third star in this binary--single encounter is not stored in the collision.log file. The only way to get information about the third star is to find this binary-single encounter in the initial.binint.log file (can be identified easily using the encounter time (here t=0.00266) and also cross-checking the id numbers for the two stars listed in the collision file).



initial.semergedisrupt.log
--------------------------

This file lists all stellar mergers that occur through binary evolution in each simulation. 

==============================  =====================================================
``t``						     time
``interaction type``		     e.g., disrupted1, disrupted2, disrupted both
``idr(mr)``						 ID_remnant(mass of the remnant)
``id1(r1)``					     ID_1 (mass of body_1)
``id2(m2)``					 	 ID_2 (mass of body_2)
``r``						 
``typer``			             merger stellar type
``type1``					     body_1 stellar type
``type2``						 body_2 stellar type
==============================  =====================================================


.. _escfile:

initial.esc.dat
---------------

As the result of dynamical encounters (and other mechanisms such as cluster tidal truncation) single stars and binaries often become unbound from the cluster potential and are ejected from the system. When this happens, the ejection is recorded in initial.esc.dat. 

==============================  =====================================================
``tcount``						     time count
``t``		     					 time
``m``						         mass
``r``					             radius
``vr``					 	         radial velocity
``vt``						 		 tangential velocity
``r_peri``			            
``r_apo``					     
``Rtidal``	
``phi_rtidal`` potential at tidal radius of the cluster 
``phi_zero``   potential at center
``E``          total energy
``J``          total ang mom
``id``
``binflag``
``m0``                           [:math:`M_{\odot}`]
``m1``                           [:math:`M_{\odot}`]
``id0``
``id1``
``a``
``e``
``startype``
``bin_startype0``	
``bin_startype1``
``rad0``
``rad1``
``tb``                          Binary orbital period [days]
``lum0``
``lum1``
``massc0``                      core mass
``massc1``
``radc0``                       core radius
``radc1``
``menv0`` 
``menv1``
``renv0``
``renv1``
``tms0``
``tms1``
``dmdt0`` mass acreting rate
``dmdt1``
``radrol0``
``radrol1``
``ospin0``  spin angular momentum of the star
``ospin1``
``B0`` mag field (Gauss)
``B1``
``formation0`` supernova (core collapsed pair instability, etc cosmic website)
``formation1``
``bacc0`` 
``bacc1``
``tacc0``
``tacc1``
``mass0_0`` initial mass 
``mass0_1``
``epoch0``  ??
``epoch1``
``bhspin`` black hole spin for the isolated bh
``bhspin1``
``bhspin2``
``ospin``
``B``     isolated mag
``formation``	isolated star			 
==============================  =====================================================


initial.morepulsars.dat
-----------------------

This files contains detailed information on all neutron stars for each simulation. For further information on treatment of neutron stars, see Ye et al. 2019, ApJ.
id 1 -100 isolated


if id 0 or id 1 binary
==============================  =====================================================
``tcount``						     time count			 
`` ``
``TotalTime ``
``binflag``
``id0 ``
``id1``
``m0``                                [:math:`M_{\odot}`]
``m1``                                [:math:`M_{\odot}`]
``B0``                                [G]
``B1``                                [G]
``P0``                                spin period[sec]
``P1``                                [sec]
``startype0``
``startype1``
``a``                                 [AU]
``ecc``
``radrol0``                           Roch ratio if >1 mass transfering
``radrol1``
``dmdt0``                            mass transfer per time
``dmdt1``
``r``                                position to cluster center
``vr``
``vt``
``bacc0``                            mass accreated to star
``bacc1``
``tacc0``
``tacc1``                            time spending acceting mass

==============================  =====================================================



initial.relaxation.dat
----------------------


==============================  =====================================================
``time``
``thetase>1.5708:f``
``q``
``<M>``
``<r>``
==============================  =====================================================

 initial.lightcollision.log
---------------------------

==============================  =====================================================
``time``
``k``
``id``
``m``
``type``
``rad``
``Eb``
``ecc``
``a``							   [AU]
``rp``                             [AU]
==============================  =====================================================

==========
Snapshots
==========
.. note::

        See :ref:`here <snapshotting>` for how to set the various snapshot parameters in the ini file
        
There are three different kinds of snapshots that CMC saves:

 * **<output>.snapshot.h5** -- every star and binary, saved every ``SNAPSHOT_DELTACOUNT`` number of code timesteps
 * **<output>.window.snapshot.h5** -- every star and binary, saved in uniform physical timesteps (set in ``SNAPSHOT_WINDOWS``)
 * **<output>.bhsnapshot.h5** -- same as <output>.shapshot.h5, but just for black holes 

Each snapshot is saved as a table in the respective hdf5 file.  To see the names of the snapshots, use ``h5ls``:

.. code-block:: bash

        h5ls <output>.window.snapshot.h5

On the window snapshots from our test example, this shows two snapshots

.. code-block:: bash

        0(t=0Gyr)                Dataset {100009/Inf}
        1(t=0.010001767Gyr)      Dataset {99233/Inf}

For the windows, this shows the number of the snapshot, and the time that the snapshot was made (in whatever units the window is using).  For the other snapshots, the time is the time in code units.

The snapshots themselves are designed to be imported as pandas tables, which each table name referring to a key in the hdf5 file.  To read in the snapshot at 10Myr:

.. ipython:: python

        import pandas as pd 
        snap = pd.read_hdf('source/example_output/output.window.snapshots.h5',key='1(t=0.010001767Gyr)')
        print(snap)

This contains all the necessary information about the state of every star and binary at this given time.  We can also see the column names

.. ipython:: python

        print(snap.columns) 

You may notice, however, that these columns are exactly the same as those in the :ref:`**<output>.esc.dat** <escfile>`  file!

The following are the columns in the snapshots but not in the escape file.
==============================  =====================================================
``luminosity``                  Luminosity of isolated stars [LSUN]
``radius``                      Radius of isolated stars [RSUN]
``bin_star_lum0``               Same as lum0
``bin_star_lum1``               Same as lum1
``bin_star_radius0``            [RSUN]
``bin_star_radius1``            [RSUN]
``bin_Eb``                      Binary binding energy
``eta``                         ?
``star.phi``                    Potential at the star's position r
==============================  =====================================================

====================
Cluster Observables
====================

cmctoolkit
__________

The `cmctoolkit <https://github.com/NicholasRui/cmctoolkit>`_ is a seperate python package specifically designed to analyze CMC snapshots.  It then computes many of the relevant astrophysical profiles of interest to observers (e.g. surface brightness profiles, number density profiles, velocity dispersions, mass-to-light ratios) allowing CMC to be directly compared to globular clusters and super star clusters in the local universe.  This is accomplished by a rigorous statistical averaging of the individual cluster orbits for each star; see `Rui et al., (2021) <https://ui.adsabs.harvard.edu/abs/2021arXiv210305033R/abstract>`_ for details.

By default, the cmctoolkit will import the last snapshot in an hdf5 snapshot file:

.. ipython:: python

        import cmctoolkit as cmct
        last_snap = cmct.Snapshot(fname='source/example_output/output.window.snapshots.h5',
                                  conv='source/example_output/output.conv.sh',
                                  dist=15, # distance to cluster in kpc
                                  z=0.0017) # metallicity 
                                 

But any snapshot in the file can be loaded by specifying the hdf5 key:

.. ipython:: python

        import cmctoolkit as cmct
        first_snap = cmct.Snapshot(fname='source/example_output/output.window.snapshots.h5',
                                  conv='source/example_output/output.conv.sh',
                                  snapshot_name='0(t=0Gyr)',
                                  dist=15, # distance to cluster in kpc
                                  z=0.0017) # metallicity
                                

As an example of what the `cmctoolkit` can do, we can create U-band surface brightness profiles for both snapshots, seeing they change due to the combined effects of stellar evolution and the early expansion of the cluster due to mass loss:

.. ipython:: python

        first_snap.add_photometry('source/output/filt_index.txt');

        u_bincenter_first, u_profile_first = first_snap.make_smoothed_brightness_profile('U', bins=80,
                                                                       min_mass=None, max_mass=None,
                                                                       max_lum=None, fluxdict=None,
                                                                       startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                                       min_logr=-1.5)

        last_snap.add_photometry('source/output/filt_index.txt');

        u_bincenter_last, u_profile_last = last_snap.make_smoothed_brightness_profile('U', bins=80,
                                                                       min_mass=None, max_mass=None,
                                                                       max_lum=None, fluxdict=None,
                                                                       startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                                       min_logr=-1.5)
        plt.plot(u_bincenter_first, u_profile_first, lw=2, label='0 Myr');
        plt.plot(u_bincenter_last, u_profile_last, lw=2, label='10 Myr');

        plt.legend(loc='lower left',fontsize=14);
        plt.xlabel('$r$ (arcsec)',fontsize=15);
        plt.ylabel('$\Sigma$ (mag/pc)',fontsize=15);
        plt.xscale('log')
        plt.xlim(5e-1, 1e3);
        @savefig plot_sbp.png width=7in
        plt.ylim(33, 5); 

See the documentation on the `cmctoolkit` for more details 