# Creates a Plummer Sphere with N=10^4 particles, using COSMIC's initial condition generators
# This imports the version of COSMIC that was downloaded and installed with CMC
from cosmic.sample import InitialCMCTable

# Generate the Singles and Binaries Pandas tables.  Note that this does not add any realistic stellar parameters
# (e.g. stellar masses, radii); these are just equal-mass point particles
Singles, Binaries = InitialCMCTable.sampler('cmc_point_mass', cluster_profile='plummer', size=10000, r_max=100)

# Scale the Cluster to Henon units (G = 1, M_cluster = 1, E_cluster = -0.25)
# Note that this option is automatically done in InitialCMCTable.write if the cluster
# isn't already normalized, but we do it here explicitly.
InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries)

# Save them to an hdf5 file for CMC
InitialCMCTable.write(Singles, Binaries, filename="plummer.hdf5")
