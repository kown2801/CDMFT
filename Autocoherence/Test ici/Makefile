CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++14
CPPFLAGS += -DNDEBUG -DHAVE_MPI
CPPFLAGS += -Iinclude 
LDFLAGS += -Llib 
LIBS += -ljson_spirit -DJSON_SPIRIT_MVALUE_ENABLED

FILES =  Patrick/Integrators.h  Patrick/Plaquette/Plaquette.h Patrick/Utilities.h 
FILES+= Patrick/Flavors.h Patrick/Hyb.h Patrick/IO.h newIO.h

all: GFULL

CDMFT: CDMFT.cpp $(FILES)
	g++ $(CPPFLAGS) $(CXXFLAGS) -o $@ CDMFT.cpp $(LDFLAGS) $(LIBS)

GFULL: GFull.cpp $(FILES)
	g++ $(CPPFLAGS) $(CXXFLAGS) -o $@ GFull.cpp $(LDFLAGS) $(LIBS)



clean: rm CDMFT GFULL


