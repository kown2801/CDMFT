CXXFLAGS += -O3 -g -fexceptions -Wall 
CPPFLAGS += -I$(HOME)/local/include/json_spirit -DNDEBUG -DHAVE_MPI
 
LDFLAGS += -L$(HOME)/local/lib -L$(HOME)/local/lib/openblas 
LIBS += -ljson_spirit -lopenblas

HEADERS_IS = Bath.h Utilities.h Hyb.h Green.h Link.h Trace.h MarkovChain.h MonteCarlo.h
HEADERS_IS+= MPIUtilities.h  

all:     IS

IS:  IS.C $(HEADERS_IS)
	mpic++ $(CPPFLAGS) $(CXXFLAGS) -o $@ IS.C $(LDFLAGS) $(LIBS)	









 
