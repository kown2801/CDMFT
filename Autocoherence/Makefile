CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++14
CPPFLAGS += -I/home/kowalski/local/include/json_spirit -DNDEBUG -DHAVE_MPI -lopenblas
CPPFLAGS += -I/home/kowalski/local/include
LDFLAGS += -L/home/kowalski/local/lib -L/home/kowalski/local/lib/openblas
LIBS += -ljson_spirit

FILES =  Patrick/Integrators.h  Patrick/Plaquette/Plaquette.h Patrick/Utilities.h 
FILES+= Patrick/Flavors.h Patrick/Hyb.h Patrick/IO.h newIO.h

all: CDMFT

CDMFT: CDMFT.cpp $(FILES)
	g++ $(CPPFLAGS) $(CXXFLAGS) -o $@ CDMFT.cpp $(LDFLAGS) $(LIBS)

GFULL: GFull.cpp $(FILES)
	g++ $(CPPFLAGS) $(CXXFLAGS) -o $@ GFull.cpp $(LDFLAGS) $(LIBS)



clean: rm CDMFT GFULL


