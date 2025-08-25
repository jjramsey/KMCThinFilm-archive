include make.inc

CPPFLAGS_COMMON = -I$(KMC_INST)/include -I$(BOOST_ROOT)/include -I$(DCMT_ROOT)/include
LDFLAGS_COMMON = -L$(KMC_INST)/lib -Wl,-rpath,$(KMC_INST)/lib

CPPFLAGS_PARALLEL = -I$(KMC_INST)/include/KMCThinFilm/parallel $(CPPFLAGS_COMMON)
CXXFLAGS_PARALLEL = $(CXXFLAGS_COMMON)
LDFLAGS_PARALLEL =  $(LDFLAGS_COMMON) -lKMCThinFilmParallel

CPPFLAGS_SERIAL = -I$(KMC_INST)/include/KMCThinFilm/serial $(CPPFLAGS_COMMON)
CXXFLAGS_SERIAL = $(CXXFLAGS_COMMON)
LDFLAGS_SERIAL =  $(LDFLAGS_COMMON) -lKMCThinFilmSerial

TARG_NAME = testFractal

all: $(TARG_NAME)Serial $(TARG_NAME)ParallelRow $(TARG_NAME)ParallelCompact $(TARG_NAME)ParallelBad

$(TARG_NAME)ParallelRow: EventsAndActions.cpp testFractal.cpp
	$(MPICXX) $(CPPFLAGS_PARALLEL) $(CXXFLAGS_PARALLEL) -c EventsAndActions.cpp
	$(MPICXX) $(CPPFLAGS_PARALLEL) $(CXXFLAGS_PARALLEL) -c testFractal.cpp
	$(MPICXX) $(LDFLAGS_PARALLEL) -o $(TARG_NAME)ParallelRow testFractal.o EventsAndActions.o

$(TARG_NAME)ParallelCompact: EventsAndActions.cpp testFractal.cpp
	$(MPICXX) $(CPPFLAGS_PARALLEL) $(CXXFLAGS_PARALLEL) -c EventsAndActions.cpp
	$(MPICXX) $(CPPFLAGS_PARALLEL) -DUSE_COMPACT_DECOMP $(CXXFLAGS_PARALLEL) -c testFractal.cpp
	$(MPICXX) $(LDFLAGS_PARALLEL) -o $(TARG_NAME)ParallelCompact testFractal.o EventsAndActions.o

$(TARG_NAME)ParallelBad: EventsAndActions.cpp testFractal.cpp
	$(MPICXX) $(CPPFLAGS_PARALLEL) $(CXXFLAGS_PARALLEL) -c EventsAndActions.cpp
	$(MPICXX) $(CPPFLAGS_PARALLEL) -DUSE_COMPACT_DECOMP -DBAD_NSTOP $(CXXFLAGS_PARALLEL) -c testFractal.cpp
	$(MPICXX) $(LDFLAGS_PARALLEL) -o $(TARG_NAME)ParallelBad testFractal.o EventsAndActions.o

$(TARG_NAME)Serial: EventsAndActions.cpp testFractal.cpp
	$(CXX) $(CPPFLAGS_SERIAL) $(CXXFLAGS_SERIAL) -c EventsAndActions.cpp
	$(CXX) $(CPPFLAGS_SERIAL) $(CXXFLAGS_SERIAL) -c testFractal.cpp
	$(CXX) $(LDFLAGS_SERIAL) -o $(TARG_NAME)Serial testFractal.o EventsAndActions.o

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(TARG_NAME)Serial $(TARG_NAME)ParallelRow $(TARG_NAME)ParallelCompact $(TARG_NAME)ParallelBad
	rm -f testdir_parallel_row/*.dat testdir_parallel_row/*.png
	rm -f testdir_parallel_compact/*.dat testdir_parallel_compact/*.png
	rm -f testdir_parallel_bad/*.dat testdir_parallel_bad/*.png
	rm -f testdir_serial/*.dat testdir_serial/*.png
