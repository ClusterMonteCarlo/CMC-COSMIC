# Creates a King profile with N=10^5 stars/binaries and a central concentration (W_0) of 6,
# using COSMIC's initial condition generators
# This imports the version of COSMIC that was downloaded and installed with CMC
from cosmic.sample import InitialCMCTable

# Generate the Singles and Binaries Pandas tables.  Here we sample stellar masses from a Kroupa 2001
# IMF, and assign 10% of the stars binary companions.  Radii are set in COSMIC.
# See (https://cosmic-popsynth.github.io/COSMIC/runpop/index.html) for more details on the options here
Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=1., primary_model='kroupa01',
                                            ecc_model='thermal', porb_model='log_uniform', qmin=-1.0,
                                            cluster_profile='king', met=0.00017, size=100000,w_0=6,
                                            params='KingProfile.ini',
                                            seed=12345,sample_porb_first=True,set_radii_with_BSE=True)

# Scale the Cluster to Henon units (G = 1, M_cluster = 1, E_cluster = -0.25)
# Note that this option is automatically done in InitialCMCTable.write if the cluster
# isn't already normalized, but we do it here explicitly.
InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries)

# Save them to an hdf5 file for CMC
InitialCMCTable.write(Singles, Binaries, filename="king.hdf5")
