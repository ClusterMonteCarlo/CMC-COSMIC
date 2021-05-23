# Creates a King profile with N=10^5 stars/binaries and a central concentration (W_0) of 6,
# using COSMIC's initial condition generators
# This imports the version of COSMIC that was downloaded and installed with CMC
# See https://clustermontecarlo.github.io/CMC-COSMIC/initialconditions/index.html for more details
from cosmic.sample import InitialCMCTable

# Generate the Singles and Binaries Pandas tables.  Here we sample stellar masses from a Kroupa 2001
# IMF, and assign 10% of the stars binary companions.  Radii are set in COSMIC.
# See (https://cosmic-popsynth.github.io/COSMIC/runpop/index.html) for more details on the options here
Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.1, primary_model='kroupa01',
                                            ecc_model='thermal', porb_model='log_uniform', qmin=0.1,
                                            cluster_profile='king', met=0.00017, zsun=0.017, size=100000,w_0=6,
                                            ,seed=12345, virial_radius=1,tidal_radius = 1e6)

# Save them to an hdf5 file for CMC with a virial radius of 1pc 
# and an initial tidal radius of 10^6 virial radii
InitialCMCTable.write(Singles, Binaries, filename="king.hdf5")
