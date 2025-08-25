# These paths should be set to the actual paths needed for compilation
KMC_INST = /data/Projects/KMC-film-growth/test-code/KMC-prototype/inst-test3
BOOST_ROOT = /data/Projects/Boost
DCMT_ROOT = /data/Projects/KMC-film-growth/ExtLibs/dcmt0.6.2

# Set the compiler and MPI wrapper to one of those available at one's workstation or cluster.
CXX = g++
MPICXX = mpicxx

# Options for the compiler
CXXFLAGS_COMMON = -Wall -O3
