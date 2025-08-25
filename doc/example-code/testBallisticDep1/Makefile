include make.inc

CPPFLAGS = -I$(KMC_INST)/include -I$(KMC_INST)/include/KMCThinFilm/serial -I$(BOOST_ROOT)/include
LDFLAGS =  -L$(KMC_INST)/lib -Wl,-rpath,$(KMC_INST)/lib -lKMCThinFilmSerial

TARG_NAME = testBallisticDep

all: $(TARG_NAME)

InitLattice.o: InitLattice.cpp InitLattice.hpp EventsAndActions.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c InitLattice.cpp

EventsAndActions.o: EventsAndActions.cpp EventsAndActions.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c EventsAndActions.cpp

$(TARG_NAME): EventsAndActions.o InitLattice.o testBallisticDep.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c testBallisticDep.cpp
	$(CXX) $(LDFLAGS) -o $(TARG_NAME) testBallisticDep.o EventsAndActions.o InitLattice.o

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(TARG_NAME)
	rm -f testdir/*.3D


