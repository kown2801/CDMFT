CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++11 
CPPINCLUDES += -I${HOME}/local/include/json_spirit -I${HOME}/local/include/openblas
CPPFLAGS += -DNDEBUG -DHAVE_MPI
LDFLAGS += -L${HOME}/local/lib
LIBS += -ljson_spirit -lopenblas -lpthread

HEADERS_IS = Bath.h Utilities.h Hyb.h Green.h Link.h Trace.h MarkovChain.h MonteCarlo.h
HEADERS_IS+= MPIUtilities.h  

all:     IS

IS:  IS.C $(HEADERS_IS)
	mpic++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ IS.C $(LDFLAGS) $(LIBS)	

clean:
	rm IS