include make.inc

CPPFLAGS = -I$(KMC_INST)/include -I$(KMC_INST)/include/KMCThinFilm/serial -I$(BOOST_ROOT)/include
LDFLAGS =  -L$(KMC_INST)/lib -Wl,-rpath,$(KMC_INST)/lib -lKMCThinFilmSerial

TARG_NAME = testFractal

all: $(TARG_NAME)

EventsAndActions.o: EventsAndActions.cpp EventsAndActions.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c EventsAndActions.cpp

$(TARG_NAME): EventsAndActions.o testFractal.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c testFractal.cpp
	$(CXX) $(LDFLAGS) -o $(TARG_NAME) testFractal.o EventsAndActions.o

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(TARG_NAME)
	rm -f testdir/*.dat testdir/*.png
