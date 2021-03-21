CXXFLAGS += -g -rdynamic -O3 -Wall -std=c++14 
CPPINCLUDES += -I${HOME}/local/include
CPPFLAGS += -DNDEBUG -DHAVE_MPI
LDFLAGS += -L${HOME}/local/lib
LIBS += -lopenblas -lpthread

HEADERS_IS = Bath.h Utilities.h Hyb.h Green.h Link.h Trace.h MarkovChain.h MonteCarlo.h IO.h
HEADERS_IS+= MPIUtilities.h  

all:     IS

IS:  IS.C $(HEADERS_IS)
	source ../scripts/export.sh > /dev/null 2>&1 ; mpic++ $(CPPINCLUDES) $(CPPFLAGS) $(CXXFLAGS) -o $@ IS.C $(LDFLAGS) $(LIBS)	

clean:
	rm IS
