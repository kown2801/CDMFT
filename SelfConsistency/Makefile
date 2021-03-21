CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++14 
CPPINCLUDES += -I${HOME}/local/include
CPPFLAGS += -DNDEBUG -DHAVE_MPI
LDFLAGS += -L${HOME}/local/lib
LIBS += -lopenblas -lpthread

FILES =  Patrick/Integrators.h  Patrick/Plaquette/Plaquette.h Patrick/Utilities.h Patrick/Hyb.h Patrick/Flavors.h
FILES+= IO.h

all: CDMFT GFULL STIFFNESS 

CDMFT: CDMFT.cpp $(FILES)
	source ../scripts/export.sh > /dev/null 2>&1; g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ CDMFT.cpp $(LDFLAGS) $(LIBS)

GFULL: GFull.cpp $(FILES)
	source ../scripts/export.sh > /dev/null 2>&1 ; g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ GFull.cpp $(LDFLAGS) $(LIBS)

STIFFNESS: STIFFNESS.cpp $(FILES)
	source ../scripts/export.sh > /dev/null 2>&1 ; g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ STIFFNESS.cpp $(LDFLAGS) $(LIBS)

clean: 
	rm CDMFT GFULL STIFFNESS
