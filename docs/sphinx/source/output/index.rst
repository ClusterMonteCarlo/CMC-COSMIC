.. _output:i

####################
Analyzing CMC Output
####################

==========
Log Files
==========
FOR FULYA



initial.conv.sh
---------------

CMC uses the standard N-body units `(He ́non 1971a) <https://link.springer.com/article/10.1007/BF00649159>`_ in all numerical calculations: the gravitational constant :math:`{G=1}`, the initial total mas :math:`{M=1}`, and the initial total energy :math:`{E_0=-1/4}`. The conversion from code units to physical units for time is done by calculating the initial half-mass relaxation time in years. The ``file initial.conv.sh`` contains all conversion factors used to convert from code to physical units. For example, to convert from code units of time to Myr, multiply by the factor timeunitsmyr, and so on. 

One can read the `initial.conv.sh` file, store the physical units and then call them as::

>>> def read_units():
>>>		s = 'initial.conv.sh'
>>>		f=open(s,'r')
>>>		count=0
>>>		while count<10:;
>>>			a=f.readline()
>>>			b=a.split('=')
>>>			if b[0]=='massunitcgs':
>>>				massunitcgs = float(b[1])
>>>				count+=1
>>>			elif b[0]=='massunitmsun':
>>>				massunitmsun = float(b[1])
>>>				count+=1
>>>			elif b[0]=='mstarunitcgs':
>>>				mstarunitcgs = float(b[1])
>>>				count+=1
>>>			elif b[0]=='mstarunitmsun':
>>>				mstarunitmsun = float(b[1])
>>>				count+=1
>>>			elif b[0]=='lengthunitcgs':
>>>				lengthunitcgs = float(b[1])
>>>				count+=1
>>>			elif b[0]=='lengthunitparsec':
>>>				lengthunitparsec = float(b[1])
>>>				count+=1
>>>			elif b[0]=='timeunitcgs':
>>>				timeunitcgs = float(b[1])
>>>				count+=1
>>>			elif b[0]=='timeunitsmyr':
>>>				timeunitsmyr = float(b[1])
>>>				count+=1
>>>			elif b[0]=='nbtimeunitcgs':
>>>				nbtimeunitcgs = float(b[1])
>>>				count+=1
>>>			elif b[0]=='nbtimeunitsmyr':
>>>				nbtimeunitsmyr = float(b[1])
>>>				count+=1
>>>		f.close()
>>>		units = []
>>>		unittype = [('m_cgs', float), ('m_msun', float), ('mstar_cgs', float), ('mstar_msun', float), ('l_cgs', float), ('l_pc', float), ('t_cgs', float),('t_myr', float), ('nbt_cgs', float), ('nbt_myr', float)]
>>>		units.append((massunitcgs, massunitmsun, mstarunitcgs, mstarunitmsun, lengthunitcgs, lengthunitparsec, timeunitcgs, timeunitsmyr, nbtimeunitcgs, nbtimeunitsmyr))
>>>		units = array(units, dtype=unittype)	
>>>		return units
>>> units=read_units() 
>>> km = units[0]['l_cgs']*1.0e-5
>>> time_units = units[0]['t_myr']
>>> m_cgs = units[0]['l_cgs']*1.0e-2
>>> kg = units[0]['m_cgs']*1.0e-3
>>> time_cgs = units[0]['t_cgs']
>>> nbt = units[0]['nbt_cgs']
>>> M_total_units = units[0]['m_msun']


initial.dyn.dat
---------------
The `initial.dyn.dat` files contains various theoretical quantities pertaining to the dynamical evolution of your cluster, such as theoretical core radius, total mass, central density, central velocity dispersion, etc. This file may be used to create, for example, plots of core radii or half-light radii versus time for a given simulation. 


================  =====================================================
``t``               time
``dt``              timestep
``tcount``
``N``
``M``
``VR``
``N_c``
``r_c``
``r_max``
``Etot``
``KE``
``PE``
``Etot_int``
``Etot_bin``
``E_cenma``
``Eesc``
``Ebesc``
``Eintesc``
``Eoops``
``Etot+Eoops``
``r_h``
``rho_0``
``rc_spitzer``
``v0_rms``
``rc_nb``
``DMse``			[:math:`{M_{\odot}}`]
``DMrejuv`` 		[:math:`{M_{\odot}}`]
``N_c_nb``
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
>>> input: type=binary m0=0.211113 m1=0.148022 R0=0.22897 R1=0.170412 Eint1=0 Eint2=0 id0=33128 id1=1255329 a=0.0908923 e=0.0641548 ktype1=0 ktype2=0
>>> status: DE/E=-1.79889e-08 DE=1.71461e-10 DL/L=2.54957e-08 DL=8.18406e-10 DE_GW/E=-0 DE_GW=0 v_esc_cluster[km/s]=77.9847 tcpu=0.01
>>> outcome: nstar=3 nobj=2:  0 [1 2] (single-binary)
>>> output: type=single m=0.0284732 R=0.215538 Eint=0 id=170307 ktype=0
>>> output: type=binary m0=0.211113 m1=0.148022 R0=0.22897 R1=0.170412 Eint1=0 Eint2=0 id0=33128 id1=1255329 a=0.09094 e=0.123848 ktype1=0 ktype2=0
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
``a``							semi-major axis [AU]
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


initial.BBHmerger.dat
---------------------

This file lists masses and merger times for all binary BH mergers occurring in the given simulation. Also listed are merger channels (ejected, in-cluster, etc.) as described in Section 9 of Kremer et al. 2019. 



initial.collision.log
---------------------

This file lists stellar types and properties for all stellar collisions occurring in the given simulation. See Sections 6 and 7 of Kremer et al. 2019 for further detail. 


initial.semergedisrupt.log
--------------------------

This file lists all stellar mergers that occur through binary evolution in each simulation. 


initial.esc.dat
---------------

As the result of dynamical encounters (and other mechanisms such as cluster tidal truncation) single stars and binaries often become unbound from the cluster potential and are ejected from the system. When this happens, the ejection is recorded in initial.esc.dat. 


initial.morepulsars.dat
-----------------------

This files contains detailed information on all pulsars for each simulation. For further information on treatment of pulsars, see Ye et al. 2019, ApJ.

==========
Snapshots
==========
FOR CARL

====================
Cluster Observables
====================

cmctoolkit
__________
FOR NICHOLAS (MAYBE)

fresca
______
FOR CARL (MAYBE)
