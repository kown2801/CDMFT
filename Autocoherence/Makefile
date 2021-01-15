CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++14 
CPPINCLUDES += -I/home/kowalski/local/include/json_spirit -I/home/kowalski/local/include/openblas
CPPFLAGS += -DNDEBUG -DHAVE_MPI
LDFLAGS += -L/home/kowalski/local/lib
LIBS += -ljson_spirit -lopenblas -lpthread

FILES =  Patrick/Integrators.h  Patrick/Plaquette/Plaquette.h Patrick/Utilities.h 
FILES+= Patrick/Flavors.h Patrick/Hyb.h Patrick/IO.h newIO.h

all: CDMFT GFULL STIFFNESS SUSCEPTIBILITY

CDMFT: CDMFT.cpp $(FILES)
	g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ CDMFT.cpp $(LDFLAGS) $(LIBS)

GFULL: GFull.cpp $(FILES)
	g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ GFull.cpp $(LDFLAGS) $(LIBS)

STIFFNESS: STIFFNESS.cpp $(FILES)
	g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ STIFFNESS.cpp $(LDFLAGS) $(LIBS)

SUSCEPTIBILITY: SUSCEPTIBILITY.cpp $(FILES)
	g++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ SUSCEPTIBILITY.cpp $(LDFLAGS) $(LIBS)
clean: 
	rm CDMFT GFULL STIFFNESS
