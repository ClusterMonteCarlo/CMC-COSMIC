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


Here we list all of the output files that are generated in a typical CMC run, 
and offer some suggestsions on how to analyze them.  Generally, the output 
comes in one of two types: `log files`, where individual events (e.g. stellar 
collisions, mergers, fewbody interactions) are logged whenenver they happen or 
when specific quantities are recorded every timestep (e.g. the number of black 
holes, the core radius, the lagrange radii), and `snapshots`, where the state 
of every star and binary in the cluster is recorded.  These are saved in plain 
text dat files and hdf5 files, respectively.

Most of the log files are in plain text (.dat) with easily readable columns 
that can be loaded in python with either ``loadtxt`` or the appropriate pandas 
commands.  However, a handful of the event log files, such as the binary 
interaction (binint.log), collision (collision.log), and stellar evolution 
distruption (semergedistrupt.log), are somewhat more difficult to parse.  We 
have provided a python package to make parsing these easier; it is located in 
the ``CMC-COSMIC/tools`` directory of the main CMC folder (the `cmc_parser.py 
<https://github.com/ClusterMonteCarlo/CMC-COSMIC/tree/master/tools>`_ file).

.. _cmcparser:

The log files can be imported using (specifying only the prefix that you provided for the output when running CMC):

.. ipython:: python

        import cmc_parser as cp
        units,binints,semergers,collisions = cp.load_ineraction_files('source/example_output/output')


We will discuss each of these below.

The HDF5 snapshots are designed to be importable as Pandas tables and are 
described below.  These snapshots can also be converted into various 
observational quantities, such as surface brightness profiles, velocity 
dispersion profiles, and more, using the `cmctoolkit` package (`Rui et al., 
2021 <https://ui.adsabs.harvard.edu/abs/2021arXiv210305033R/abstract>`_).  We 
include the python code from the main package in the ``CMC-COSMIC/tools`` 
folder, but please see `here <https://github.com/NicholasRui/cmctoolkit>`_ for 
complete instructions and useage.

.. note::

        Here we are using the output from the :ref:`realistic-initial-conditions` example, running it for 15 Myr (though note that the default ini file runs it for 13.8 Gyr) and having it print a window snapshot every 10 Myr.  On 28 cores of the Vera machine, this takes about 10 minutes.  

==========
Code Units
==========

CMC uses the standard N-body units `(Hénon 1971a) 
<https://link.springer.com/article/10.1007/BF00649159>`_ in all numerical 
calculations: the gravitational constant :math:`{G=1}`, the initial total mass 
:math:`{M=1}`, and the initial total energy :math:`{E_0=-1/4}`. The conversion 
from code units to physical units for time is done by calculating the initial 
half-mass relaxation time in years. The file ``<output>.conv.sh`` contains all 
conversion factors used to convert from code to physical units. For example, to 
convert from code units of time to Myr, multiply by the factor timeunitsmyr, 
and so on. **Unless otherwise noted, all values below are in code units.** 

This file can be parsed using the above command, or using a stand-alone command for only the conversion dictionary:

.. ipython:: python

        import cmc_parser as cp
        conv_file = cp.conversion_file('source/example_output/output.conv.sh')
        print(conv_file.mass_msun)
        print(conv_file.__dict__)

==========
Log Files
==========

Listed below are the output files of CMC with the variables printed out. 

initial.dyn.dat
---------------

The **initial.dyn.dat** files contains various theoretical quantities pertaining 
to the dynamical evolution of your cluster, such as theoretical core radius, 
total mass, central density, central velocity dispersion, etc. This file may be 
used to create, for example, plots of core radii or half-mass radii versus time 
for a given simulation. Be warned that these theoretical quantities are not 
necessarily the same as the way observers would define these quantities. 
However, these quantities are all useful for establishing a general idea of how 
the cluster evolves.

================  =====================================================
``t``               Time [code units]
``dt``              Time step
``tcount``          Time count
``N``               Total number of objects (single+binary)
``M``               Total mass (in units of initial mass)
``VR``              Measure how far away the cluster is from virial equilibrium (-2.0*Etotal.K/Etotal.P)
``N_c``             Total number of stars within core radius
``r_c``             Core radius
``r_max``           Maximum radius of a star 
``Etot``            Total energy 
``KE``              Total kinetic energy 
``PE``              Total potential energy 
``Etot_int``        Total internal energy of single stars
``Etot_bin``        Total internal energy of binary stars
``E_cenma``         Central BH mass (i.e., when there is a central IMBH)
``Eesc``            Total energy of the escaped single stars
``Ebesc``           Total energy of the escaped binary stars
``Eintesc``         Total internal energy in the escaped stars
``Eoops``           Energy error loss due to Stodolkiwecz's potential correction 
``Etot+Eoops``      Total energy + Eoops
``r_h``             Half-mass radius
``rho_0``           Core density
``rc_spitzer``      Core radius as defined in Spitzer 1987: :math:`\sqrt{3  \sigma_0^2}{4 \pi \rho_0}`
``v0_rms``          Rms velocity dispersion at the cluster center
``rc_nb``           Core radius calculated with density weighted averages as in Casertano & Hut (1985)
``DMse``            Total mass loss from the cluster per time step due to stellar evolution [:math:`{M_{\odot}}`]
``DMrejuv`` 	     Mass loss from rejuvenation per time step [:math:`{M_{\odot}}`]
``N_c_nb``          Number of stars within the core: :math:`\frac{4 \pi}{3} rc_{\rm nb}^3  \frac{n_{\rm c}}{2}`
================  =====================================================

initial.binint.log
------------------

Over the course of the evolution of the cluster, single stars and binaries will 
frequently undergo three- and four-body dynamical encounters, which are 
integrated directly in CMC using the Fewbody package (Fregeau et al. 2007). The 
file **initial.binint.log** records all input and output parameters (e.g., 
component masses, IDS, stellar types, semi-major axes, etc.) each time fewbody 
is called. 

Every encounter information is printed between two lines of asterisks.
Below is an exemplary output:

.. code-block:: bash

      ********************************************************************************
      type=BS t=5.85010072e-06
      params: b=1.46611 v=0.379587
      input: type=single m=0.0284732 R=0.215538 Eint=0 id=170307 ktype=0
      input: type=binary m0=0.211113 m1=0.148022 R0=0.22897 R1=0.170412 Eint1=0 Eint2=0 id0=33128 id1=1255329 a=0.0908923 e=0.0641548 ktype1=0 ktype2=0 status:      DE/E=-1.79889e-08 DE=1.71461e-10 DL/L=2.54957e-08 DL=8.18406e-10 DE_GW/E=-0 DE_GW=0 v_esc_cluster[km/s]=77.9847 tcpu=0.01
      outcome: nstar=3 nobj=2:  0 [1 2] (single-binary)
      output: type=single m=0.0284732 R=0.215538 Eint=0 id=170307 ktype=0
      output: type=binary m0=0.211113 m1=0.148022 R0=0.22897 R1=0.170412 Eint1=0 Eint2=0 id0=33128 id1=1255329 a=0.09094 e=0.123848 ktype1=0 ktype2=0
      ********************************************************************************

==============================  =====================================================
``type``						         Encounter type (BS for binary-single or BB for binary-binary)
``t``							         Encounter time
``b``							         Impact parameter [units of :math:`a` for binary-single or :math:`a_1+a_2` for binary-binary]
``v``							         Relative velocity at infinity [:math:`v_c`]
``m``							         Mass [:math:`{M_{\odot}}`]
``R``							         Radius [:math:`R_{\odot}`]
``Eint``			                  Internal energy
``id``						         ID number 
``kytpe``					         Stellar type
``a``							         Semi-major axis [AU]
``e``							         Eccentricity
``dE/E``			                  Fractional change in energy
``DE``                           Total change in energy
``DL/L``                         Fractional change in angular momentum
``DL``                           Change in angular momentum
``DE_GW``                        Energy loss due gravitational wave emission
``v_esc_cluster``			         Escape speed of the cluster where the encounter occured [km/s]
``tcpu``                         CPU time for integration (usually ~milliseconds, unless it's a GW capture)
``nstar``					         Number of stars
``nobj``						         Number of objects (single/binary)
``i [j k]``					         Final configuration after encounter, e.g.,  0 [1 2] (single-binary)
==============================  =====================================================

Objects are labelled starting from 0 to 3. The binary-single and binary-binary 
encounters are denoted as BS and BB, respectively. For type=binary, indices 0 
and 1 in mass, radius,id,etc. denote the primary and secondary objects in a 
binary.

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

0:1 denotes the fact that objects 0 and 1 have merged, and [0 1] indicates that 
objects 0 and 1 have formed a binary. The same is true for any pairs from 0 to 
3.

While the `binint` file is easy to read, it can be difficult to parse.  Using 
the ``load_interaction_files`` command from :ref:`above <cmcparser>` provides 
the ``binints`` object, a python list of dictionaries of every encounter:

.. ipython:: python

        print(binints[0].__dict__)

        # input binaries is a list that can be printed with:
        print(binints[0].in_binaries[0].__dict__)

        # and the individual stars of that binary can be accessed with:
        print(binints[0].in_binaries[0].star1.__dict__)

        # so if I wanted, for instance, the radius of star two of the input binary in the first encounter:
        print(binints[0].in_binaries[0].star2.r_RSUN)

        # If I wanted to know the escape speed of the cluster where this encouner occured, I can access that with
        print(binints[0].vesc_KMS)

initial.bh.dat
-------------

This file contains the number of BHs (as well as BH binaries, etc.) at each 
dynamical time step. This is useful to plot, for example, the number of 
retained BHs versus time. For BH mergers, you want to look in 
**initial.bhmerger.dat**, which records all BH-BH mergers that occur inside the 
cluster during the cluster evolution.

==============================  =====================================================
``tcount``						      Time count
``Totaltime``					      Total time
``Nbh,tot``						      Total number of BHs
``Nbh,single``					      Number of single BHs
``Nbinarybh``					      Number of binary BHs
``Nbh-bh``						      Number of BH-BH binaries
``Nbh-nonbh``			            Number of BH-non BH binaries
``Nbh-ns``					         Number of BH-NS binaries
``Nbh-wd``						      Number of BH-WD binaries
``N_bh-star``				 	      Number of stars including MS stars and giants 
``Nbh-ms``						      Number of BH-MS binaries	
``Nbh-postms``			            Number of BH-giant binaries
``fb_bh``						      Number of binaries containing a black hole / total number of systems containing a black hole 
==============================  =====================================================

initial.bh.esc.dat
----------------
This file contains the number of ejected BHs at each dynamical time step. It 
includes the same columns in the **initial.bh.dat** file.

initial.bhmerger.dat
---------------------

List of all binary black hole mergers that occur in the cluster (note this does 
not include BBHs that may be ejected from the cluster and merge later).  There 
are four categories of mergers that occur inside the cluster:

 * **isolat-binary** - merger that occurs in a binary, but not due to GW capture
 * **binary-single** - merger that occurs due to GW capture during a binary-single encounter
 * **binary-binary** - merger that occurs due to GW capture during a binary-binary encounter
 * **single-single** - merger that occurs due to GW capture between two isoalted black holes

==============================  =====================================================
``time``                        Time merger occurs
``type``                        What kind of merger was this 
``r``                           Radius in cluster where merger occured
``id1``                         ID of primary
``id2``                         ID of secondary 
``m1``                          Mass of primary :math:`[M_{\odot}]`
``m2``                          Mass of secondary :math:`[M_{\odot}]`
``spin1``                       Spin of primary 
``spin2``                       Spin of secondary 
``m_final``                     Mass of merger remnant :math:`[M_{\odot}]`
``spin_final``                  Spin of merger remnant
``vkick``                       Kick merger remnant recieves [km/s]
``v_esc``                       Escape speed of cluster where merger occurs [km/s]
``a_final``                     Last semi-major axis recorded for binary (see note) [AU]
``e_final``                     Last eccentricity recorded for binary
``a_50M``                       (Newtonian) semi-major axis when the BHs were 50M apart; only for binary-single or binary-binary [AU] 
``e_50M``                       (Newtonian) eccentricity when BHs were 50M apart  
``a_100M``                      Same, but 100M apart [AU] 
``e_100M``                      "
``a_500M``                      "
``e_500M``                      "
==============================  =====================================================

 .. DANGER::

        The ``a_final`` and ``e_final`` parameters change depending on the type of encounter.  For binary-single and binary-binary GW captures, these record the 
        (Newtonian) semi-major axis and eccentricity at 10M (when we consider the BHs to have mergred.  However, this is an unreliable quantity, since the orbit 
        is decidedly non-Newtonian at that point.  If you want eccentricities, use ``a_100M`` and ``e_100M``, or preferably the outermost value above).

        For single-single GW captures, ``a_final`` and ``e_final`` are the semi-major axis and eccentricity that the GW capture formed at.  For isolat-binary 
        mergers, it's the last semi-major axis and eccentricity that were recorded in the cluster.


initial.collision.log
---------------------

This file lists stellar types and properties for all stellar collisions 
occurring in the given simulation. See Sections 6 and 7 of Kremer et al. 2019 
for further detail. 

==============================  =====================================================
``t``						           Collision time
``interaction type``		        Interaction type e.g., single-binary, binary-binary, etc.
``idm(mm)``						     ID_merger(mass of merged body)
``id1(m1)``					        ID_1 (mass of collided body_1)
``id2(m2)``					 	     ID_2 (mass of the collided body_2)
``r``						           Distance from the center of cluster
``typem``			              Merger stellar type
``type1``					        Stellar type of body_1
``type2``						     Stellar type of body_2 
``b``                            Impact parameter at infinity [:math:`R_{\odot}`]
``vinf``                         Relative velocity of two objects at infinity [km/s] 
``rad1``                         Radius of body_1
``rad2``                         Radius of body_2
``rperi``                        Pericenter distance at collision
``coll_mult``                    Collison multiplyer e.g., sticky sphere (``coll_mul`` = 1), TDE (``coll_mul``> 1)
==============================  =====================================================


The single-single, binary-single, etc indicate whether the collision occurred 
during a binary encounter or not. When there are three stars listed for the 
collision it means that all three stars collided during the encounter. This is 
rare, but it does happen occasionally. Typically, one will see something like:

.. code-block:: bash 

      t=0.00266079 binary-single idm=717258(mm=1.0954) id1=286760(m1=0.669391):id2=415309 (m2=0.426012) (r=0.370419) typem=1 type1=0 type2=0

In this case the colliding stars are m1=0.66 and m2=0.42. The information about 
the third star in this binary--single encounter is not stored in the 
collision.log file. The only way to get information about the third star is to 
find this binary-single encounter in the **initial.binint.log** file (can be 
identified easily using the encounter time (here t=0.00266) and also 
cross-checking the id numbers for the two stars listed in the collision file).

Similarly to the binint file, the collision file can be processed using the :ref:`load_interaction_file <cmcparser>` command

.. ipython:: python

        print(collisions[0].__dict__)


initial.semergedisrupt.log
--------------------------

This file lists all stellar mergers that occur through binary evolution in each simulation. 

==============================  =====================================================
``t``						            Time
``interaction type``		         Interaction type e.g., disrupted1, disrupted2, disrupted both
``idr(mr)``						      ID_remnant(mass of the remnant)
``id1(r1)``					         ID_1 (mass of body_1)
``id2(m2)``					 	      ID_2 (mass of body_2)
``r``						            Distance from the center of cluster
``typer``			               Stellar type of merger
``type1``					         Stellar type of body_1 
``type2``						      Stellar type of body_2 
==============================  =====================================================


The semergedisrupt file can also be processed using the :ref:`load_interaction_file <cmcparser>` command

.. ipython:: python

        print(semergers[0].__dict__)

.. _escfile:

initial.esc.dat
---------------

As the result of dynamical encounters (and other mechanisms such as cluster 
tidal truncation) single stars and binaries often become unbound from the 
cluster potential and are ejected from the system. When this happens, the 
ejection is recorded in **initial.esc.dat**. In particular, this ejection 
process plays an intimate role in the formation of merging BH binaries. If a 
BH-BH binary is ejected from the cluster with sufficiently small orbital 
separation it may merge within a Hubble time and be a possible LIGO source. To 
determine the number of such mergers, calculate the inspiral times for all 
BH-BH binaries that appear in the **initial.esc.dat** file.


Parameters with a `_0` (i.e., mass, radius, star type, etc) correspond to the 
primary star in a binary. There is also the same column for the secondary star 
with `_0` replaced by `_1` in the **initial.esc.dat** file. Parameters without 
indicies indicate single stars.  

==============================  =====================================================
``tcount``						     Time count
``t``		     					     Time
``m``						           Mass [:math:`M_{\odot}`]. If the object is binary,  ``m`` corresponds to total mass of the primary and secondary stars 
``r``					              Radius
``vr``					 	        Radial velocity
``vt``						 		  Tangential velocity
``r_peri``			              Pericenter of star's orbit in the cluster when it was ejected    
``r_apo``                       Apocenter of star's orbit in the cluster 
``Rtidal``	                    Tidal radius
``phi_rtidal``                  Potential at the tidal radius
``phi_zero``                    Potential at center
``E``                           Total energy
``J``                           Total angular momentum
``id``                          Single ID number
``binflag``                     Binary flag. If ``binflag`` = 1, the object is binary; otherwise single
``m0``                          Primary mass [:math:`M_{\odot}`]
``id0``                         Primary ID number
``a``                           Semi-major axis [AU]
``e``                           Eccentricity
``startype``                    Single star type
``bin_startype0``	              Primary star type 
``rad0``                        Primary radius [:math:`R_{\odot}`]
``tb``                          Binary orbital period [days]
``lum0``                        Primary luminosity [:math:`L_{\odot}`]
``massc0``                      Primary core mass [:math:`M_{\odot}`
``radc0``                       Primary core radius [:math:`R_{\odot}`]
``menv0``                       Primary envelope mass [:math:`M_{\odot}`]
``renv0``                       Primary envelope radius [:math:`R_{\odot}`]
``tms0``                        Primary timescale of the main sequence
``dmdt0``                       Primary mass accreting rate 
``radrol0``                     Ratio of Roche Lobe to radius
``ospin0``                      Primary spin angular momentum
``B0``                          Primary magnetic field [G]
``formation0``                  Primary formation channel for supernova, e.g., core collapse, pair instability, etc.)
``bacc0``                       Mass accreted to the primary
``tacc0``                       Time spent accreting mass to the primary 
``mass0_0``                     Primary initial mass 
``epoch0``                      
``bhspin``                      BH spin (if single)
``bhspin1``                     BH spin for primary (if binary)  
``ospin``                       Single star spin angular momentum
``B``                           Single star magnetic field [G]
``formation``	                 Single star formation channel for supernova			 
==============================  =====================================================


initial.morepulsars.dat
-----------------------

This files contains detailed information on all neutron stars for each 
simulation. For further information on treatment of neutron stars, see Ye et 
al. 2019, ApJ.

==============================  =====================================================
``tcount``						           Time count			 
``TotalTime``                         Total time
``binflag``                           Binary flag 
``id0``                               ID number
``m0``                                Mass [:math:`M_{\odot}`]
``B0``                                Magnetic field [G]
``P0``                                Spin period [sec]
``startype0``                         Star type
``a``                                 Semi-major axis[AU]
``ecc``                               Eccentricity
``radrol0``                           Roche ratio (if > 1, mass transfering)
``dmdt0``                             Mass transfer rate 
``r``                                 Distance from the cluster center
``vr``                                Radial velocity
``vt``                                Tangential velocity
``bacc0``                             Mass accreted to star
``tacc0``                             Time spent accreting mass 
==============================  =====================================================



initial.log
------------

Each time step, cluster information is printed between two lines of asterisks.
Below is an exemplary output:

.. code-block:: bash

      ******************************************************************************
      tcount=1 TotalTime=0.0000000000000000e+00 Dt=0.0000000000000000e+00
      Etotal=-0.514537 max_r=0 N_bound=1221415 Rtidal=111.234
      Mtotal=1 Etotal.P=-0.499519 Etotal.K=0.249522 VRatio=0.99905
      TidalMassLoss=0
      core_radius=0.361719 rho_core=7.18029 v_core=0.832785 Trc=994.138 conc_param=0 N_core=135329
      trh=0.100838 rh=0.811266 rh_single=0.811936 rh_binary=0.801647
      N_b=38407 M_b=0.0752504 E_b=0.26454
      ******************************************************************************
      
==============================  =====================================================
``tcount``                       Time count
``TotalTime``                    Total time
``Dt``                           Time step
``Etotal``                       Total energy
``max_r``                        Maximum radius of a star 
``N_bound``                      Number of objects bound to the cluster
``Rtidal``                       Tidal radius of the cluster
``Mtotal``                       Total mass of the cluster
``Etotal.P``                     Total potential energy of the cluster
``Etotal.K``                     Total kinetic energy of the cluster
``VRatio``                       Virial ratio
``TidalMassLoss``                Mass lost through the tidal radius
``core_radius``                  Core radius 
``rho_core``                     Core density
``v_core``                       Velocity dispersion in the core 
``Trc``                          Core relaxation timescale
``conc_param``                   King concentration parameter
``N_core``                       Number of objects within core radius
``trh``                          Half-mass relaxation time
``rh``                           Half-mass radius
``rh_single``                    Half-mass radius of single objects
``rh_binary``                    Half-mass radius of binaries
``N_b``                          Total number of binaries
``M_b``                          Total mass of binaries
``E_b``                          Total energy of binaries
==============================  =====================================================
      
      

Note that this is also printed to ``stdout`` every timestep.


initial.tidalcapture.log
------------------------

 This files contains information on tidal capture events for each simulation. 
 
* **time**                                                    - tidal capture time
* **interaction_type**                                        - (SS_COLL_GW)
* **(id1,m1,k1)+(id2,m2,k2)->[(id1,m1,k1)-a,e-(id2,m2,k2)]**  - (ID, mass and star type of interacting stars) -> [(ID, mass, stary type of the primary) - semi-major axis, eccentricity - (ID, mass, stary type of the secondary)]
 
initial.triple.dat
------------------

List of triples formed dynamically in the cluster as a result of three- and four-body dynamical encounters. 

==============================  =====================================================
``time``                         Time
``min0``                         Mass of inner object `_0` [:math:`M_{\odot}`]
``min1``                         Mass of inner object `_1` [:math:`M_{\odot}`]
``mout``                         Mass of outer object [:math:`M_{\odot}`]
``Rin0``                         Radius of inner object `_0` [:math:`R_{\odot}`]
``Rin1``                         Radius of inner object `_1` [:math:`R_{\odot}`]
``Rout``                         Radius of outer object [:math:`R_{\odot}`]
``ain``                          Semi-major axis of inner binary [AU]
``aout``                         Semi-major axis of outer binary [AU]
``ein``                          Eccentricity of inner binary
``eout``                         Eccentricity of outer binary
``ktypein0``                     Inner object `_0` stellar type
``ktypein1``                     Inner object `_1` stellar type
``kytpeout``                     Outer object stellar type
``Tlk_quad``                     Quadrupole Kozai-Lidov timescale [yr]
``Tlk_oct``                      Octupole Kozai-Lidov timescale: tlkquad/epsoct [yr]
``eps_oct``                      Octupole parameter
``T_GR``                         1PN precession timescale [yr]
``eps_GR``                       GR parameter: Tlk_quad/T_GR
==============================  =====================================================

initial.lagrad.dat
-------------------

This file contains the lagrange radii enclosing a given percentage of the cluster's 
total mass. So for example, the 10% lagrange radii printed in the 
**initial.lagrad.dat** file is the radius at a given time that encloses 10% of 
the mass. The different columns in that file give 0.1%, 5%, 99%, etc. lagrange 
radii.

initial.v2_rad_lagrad.dat
-------------------------

List of the sum of radial velocity :math:`v_{r}` within Lagrange 
radii enclosing a given percentage of the cluster's total mass.

initial.v2_tan_lagrad.dat
-------------------------

List of the sum of tangential velocity :math:`v_{t}` within 
Lagrange radii enclosing a given percentage of the cluster's total mass.


initial.nostar_lagrad.dat
------------------------

List of the number of stars within Lagrange radii enclosing a given 
percentage of the cluster's total mass.

initial.rho_lagrad.dat
---------------------

List of the density within Lagrange radii enclosing a given 
percentage of the cluster's total mass.

initial.avemass_lagrad.dat
--------------------------

List of the average mass :math:`\langle m \rangle` within Lagrange radii 
enclosing a given percentage of the cluster's total mass in units of solar mass 
[:math:`M_{\odot}`].

initial.ke_rad_lagrad.dat
------------------------

List of the total radial kinetic energy :math:`T_{r}` within 
Lagrange radii enclosing a given percentage of the cluster's total mass in code 
units.

initial.ke_tan_lagrad.dat
------------------------

List of the total tangenial kinetic energy :math:`T_{t}` within 
Lagrange radii enclosing a given percentage of the cluster's total mass in code 
units.

initial.lagrad0-0.1-1.dat
-------------------------

List of the lagrange radii for the masses in range 0.1 :math:`M_{\odot}` < m < 1 :math:`M_{\odot}`.

initial.lagrad1-1-10.dat
------------------------

List of the lagrange radii for the masses in range 1 :math:`M_{\odot}` < m < 10 :math:`M_{\odot}`.

-------------------------

initial.lagrad2-10-100.dat
--------------------------

List of the lagrange radii for the masses in range 10 :math:`M_{\odot}` < m < 100 :math:`M_{\odot}`.


initial.lagrad3-100-1000.dat
----------------------------

List of the lagrange radii for the masses in range 100 :math:`M_{\odot}` < m < 10000 :math:`M_{\odot}`.


initial.lagrad_10_info.dat
--------------------------

This file containts dynamical information of the cluster at 10 lagrange radius.

initial.core.dat
----------------

Information for the core that contains no remnants.

==============================  =====================================================
``time``                         Time
``rho_norem``                    Density of the core
``v_rms_norem``                  Velocity dispersion of the core
``rc_norem``                     Core radius
``r_spitzer_norem``              Spitzer radius 
``m_ave_norem``                  Average mass within this core radius
``n_norem``                      
``N_norem``
``T_rc_norem``
==============================  =====================================================

initial.bin.dat
--------------

This file contains information on binaries.

==============================  =====================================================
``t``                           Total time
``N_b``                         Number of binaries
``M_b``                         Mass of binaries
``E_b``                         Total energy of binaries
``r_h,s``                       Half-mass radius of single objects
``r_h,b``                       Half-mass radius of binaries
``rho_c,s``                     Core density for single objects
``rho_c,b``                     Core density for binaries
``N_bb``                        Number of binary-binary interactions
``N_bs``                        Number of binary-single interactions
``f_b,c``                       Binary fraction in the core
``f_b``                         Binary fraction
``E_bb``                         
``E_bs``
``DE_bb``
``DE_bs``
``N_bc,nb``                      Number of existing bound binaries in the core
``f_b,c,nb``                     Fraction of existing bound binaries in the core
``N_bc``                         Number of all binaries including the escaped and destroyed ones in the core
==============================  =====================================================

initial.bhformation.dat
-----------------------

This file contains information about newly formed BHs.

==============================  =====================================================
``time``                        Time of BH formation
``r``                           Position in cluster
``binary?``                     Whether binary or not at the time of stellar collapse
``ID``                          ID of BH
``zams_m``                      Mass of progenitor at t=0
``m_progenitor``                Mass of progenitor before explosion
``bh mass``                     Mass of BH
``bh_spin``                     Spin of BH
``birth-kick``                  Birth kick magnitude of natal kick [km/s]
``vsarray``                     Array of natal kicks 
==============================  =====================================================


initial.3bb.log
---------------

This file contains information of three body binaries (triples).

==============================  =====================================================
``time``                         Time
``k1``                           Stellar type of object `_1`
``k2``                           Stellar type of object `_2`
``k3``                           Stellar type of object `_3`
``id1``                          ID of object `_1`
``id2``                          ID of object `_2`
``id3``                          ID of object `_3`
``m1``                           Mass of object `_1`
``m2``                           Mass of object `_2`
``m3``                           Mass of object `_3`
``ave_local_mass``               Average local mass
``n_local``                      Local number density
``sigma_local``                  The local 3D velocity dispersion
``eta``                          Hardness of the inner binary
``Eb``                           Binding energy of the inner binary
``ecc``                          Eccentricty
``a``                            Semi-major axis [AU]
``r_peri``                       Pericenter distance [AU]
``r(bin)``                       Radial distance of the binary
``r(single)``                    Radial distance of the single object
``vr(bin)``                      Radial velocity of binary
``vt(bin)``                      Tangential velocity of binary
``vr(single)``                   Radial velocity of single object
``vt(single)``                   Tangential velocity of single object
``phi(bin)``                     Potential energy of the binary
``phi(single)``                  Potential energy of the single object
``delta_PE``                     Change of the potential energy 
``delta_KE``                     Change of the kinetic energy 
``delta_E(interaction)``         Change of the total energy per interaction
``delta_E(cumulative)``          Change of the total energy for all 3-body interactions         
``N_3bb``                        The number of triples formed
==============================  =====================================================


initial.3bbprobability.log  
--------------------------
Average rate and probability of three-body binary formation in the timestep calculated from the innermost 300 triplets of single stars considered for three-body binary formation.

==============================  =====================================================
``time``                         Time
``dt``                           Time step
``dt*N/log(gamma*N)``            
``Rate_3bb``
``P_3bb``                        Probability of binary formation
``r``                             
==============================  =====================================================


initial.lightcollision.log
---------------------------

List of stars that fail to form triples.

==============================  =====================================================
``time``                            Time
``k1``                              Variable to store most massive star index
``k2``                              Second most massive star index
``k3``                              Least massive star index
``id``                              ID number
``m``                               Mass
``type``                            Stellar type
``rad``                             Stellar radius
``Eb``                              Binding energy
``ecc``                             Eccentricity
``a``			                        Semi-major axis [AU]
``rp``                              m1 * m2 * madhoc / (eta * sqr(sigma_local)) [AU]
==============================  =====================================================


==========
Snapshots
==========
.. note::

        See :ref:`here <snapshotting>` for how to set the various snapshot parameters in the ini file
        
There are three different kinds of snapshots that CMC saves:

 * **output.snapshot.h5** -- every star and binary, saved every ``SNAPSHOT_DELTACOUNT`` number of code timesteps
 * **output.window.snapshot.h5** -- every star and binary, saved in uniform physical timesteps (set in ``SNAPSHOT_WINDOWS``)
 * **output.bhsnapshot.h5** -- same as output.shapshot.h5, but just for black holes 

Each snapshot is saved as a table in the respective hdf5 file.  To see the 
names of the snapshots, use ``h5ls``:

.. code-block:: bash

        h5ls output.window.snapshot.h5

On the window snapshots from our test example, this shows two snapshots

.. code-block:: bash

        0(t=0Gyr)                Dataset {100009/Inf}
        1(t=0.010001767Gyr)      Dataset {99233/Inf}

For the windows, this shows the number of the snapshot, and the time that the 
snapshot was made (in whatever units the window is using).  For the other 
snapshots, the time is the time in code units.

The snapshots themselves are designed to be imported as pandas tables, which 
each table name referring to a key in the hdf5 file.  To read in the snapshot 
at 10Myr:

.. ipython:: python

        import pandas as pd 
        snap = pd.read_hdf('source/example_output/output.window.snapshots.h5',key='1(t=0.010001767Gyr)')
        print(snap)

This contains all the necessary information about the state of every star and 
binary at this given time.  We can also see the column names

.. ipython:: python

        print(snap.columns) 

You may notice, however, that these columns are exactly the same as those in 
the :ref:`output.esc.dat <escfile>`  file!

The following are the columns in the snapshots but not in the escape file.

==============================  =====================================================
``luminosity``                  Luminosity of isolated stars [LSUN]
``radius``                      Radius of isolated stars [RSUN]
``bin_star_lum0``               Same as lum0
``bin_star_lum1``               Same as lum1
``bin_star_radius0``            [RSUN]
``bin_star_radius1``            [RSUN]
``bin_Eb``                      Binary binding energy
``eta``                         Binary hardness
``star.phi``                    Potential at the star's position r
==============================  =====================================================

====================
Cluster Observables
====================

The `cmctoolkit <https://github.com/NicholasRui/cmctoolkit>`_ is a seperate 
python package specifically designed to analyze CMC snapshots.  It then 
computes many of the relevant astrophysical profiles of interest to observers 
(e.g. surface brightness profiles, number density profiles, velocity 
dispersions, mass-to-light ratios) allowing CMC to be directly compared to 
globular clusters and super star clusters in the local universe.  This is 
accomplished by a rigorous statistical averaging of the individual cluster 
orbits for each star; see `Rui et al., (2021) 
<https://ui.adsabs.harvard.edu/abs/2021arXiv210305033R/abstract>`_ for details.

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
                                

As an example of what the `cmctoolkit` can do, we can create V-band surface 
brightness profiles for both snapshots, seeing they change due to the combined 
effects of stellar evolution and the early expansion of the cluster due to mass 
loss:

.. ipython:: python

        first_snap.add_photometry('source/output/filt_index.txt');

        v_bincenter_first, v_profile_first = first_snap.make_smoothed_brightness_profile('V', bins=80,
                                                                       min_mass=None, max_mass=None,
                                                                       max_lum=None, fluxdict=None,
                                                                       startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                                       min_logr=-1.5)

        last_snap.add_photometry('source/output/filt_index.txt');

        v_bincenter_last, v_profile_last = last_snap.make_smoothed_brightness_profile('V', bins=80,
                                                                       min_mass=None, max_mass=None,
                                                                       max_lum=None, fluxdict=None,
                                                                       startypes=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                                                       min_logr=-1.5)
        plt.plot(v_bincenter_first, v_profile_first, lw=2, label='0 Myr');
        plt.plot(v_bincenter_last, v_profile_last, lw=2, label='10 Myr');

        plt.legend(loc='lower left',fontsize=14);
        plt.xlabel('$r$ (arcsec)',fontsize=15);
        plt.ylabel('$\Sigma_V$ (mag/arcsec$^2$)',fontsize=15);
        plt.xscale('log')
        plt.xlim(5e-1, 1e3);
        @savefig plot_sbp.png width=7in
        plt.ylim(33, 5); 

Examples of other properties which can be easily computed by `cmctoolkit`:

 * Stellar magnitudes and colors (in blackbody limit)
 * Velocity dispersion profiles
 * Mass functions and mass function slopes
 * Binary fractions
 * Blue stragglers

See the documentation on the `cmctoolkit` for more details.
