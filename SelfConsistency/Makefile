CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++14 
CPPINCLUDES += 
CPPFLAGS += -DNDEBUG -DHAVE_MPI
LDFLAGS += -L${HOME}/local/lib
LIBS += -lopenblas -lpthread

FILES =  Patrick/Integrators.h  Patrick/Plaquette/Plaquette.h Patrick/Utilities.h Patrick/Hyb.h Patrick/Flavors.h
FILES+= IO.h nlohmann_json.hpp

all: CDMFT GFULL STIFFNESS 

CDMFT: CDMFT.cpp $(FILES)
	source ../scripts/export.sh > /dev/null 2>&1; g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ CDMFT.cpp $(LDFLAGS) $(LIBS)

GFULL: GFULL.cpp $(FILES)
	source ../scripts/export.sh > /dev/null 2>&1 ; g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ GFULL.cpp $(LDFLAGS) $(LIBS)

STIFFNESS: STIFFNESS.cpp $(FILES)
	source ../scripts/export.sh > /dev/null 2>&1 ; g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ STIFFNESS.cpp $(LDFLAGS) $(LIBS)

clean: 
	rm CDMFT GFULL STIFFNESS
