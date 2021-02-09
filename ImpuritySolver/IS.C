#include "MonteCarlo.h"
#include "IO.h"

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
#endif
	try {
		
		if(argc != 4) throw std::runtime_error("Usage : IS inputFolder outputFolder filename");
		
		std::string inputFolder = argv[1];
		std::string outputFolder = argv[2];
		std::string fileName = argv[3];

		std::time_t time;
		mpi::cout = mpi::every;
		mpi::cout << "Start task at " << std::ctime(&(time = std::time(NULL))) << std::flush;
		
		Ut::Simulation simulation(inputFolder, outputFolder, fileName);
		json_spirit::mObject const& jParams = simulation.params();
		Ma::MarkovChain* markovChain = 0;
		
		{	
			//Reading all the input files and keeping in memory the output filename

			json_spirit::mValue jHyb;
			mpi::read_json(inputFolder + jParams.at("HYB").get_str(), jHyb);
			
			json_spirit::mArray jLink;
			IO::readLinkFromParams(jLink, inputFolder,jParams);

			markovChain = new Ma::MarkovChain(jParams, jHyb.get_obj(), jLink, simulation);
		}

		MC::MonteCarlo(*markovChain, simulation);
		
		delete markovChain; markovChain = 0;
		
		mpi::cout = mpi::every;
		mpi::cout << "Task of worker finished at " << std::ctime(&(time = std::time(NULL))) << std::flush;
	}
	catch (std::exception& exc) {
		std::cerr << exc.what() << "( Thrown from worker " << mpi::rank() << " )" << std::endl;
		
#ifdef HAVE_MPI
		MPI_Abort(MPI_COMM_WORLD, -1);
#endif 
		return -1;
	}
	catch (...) {
		std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
		
#ifdef HAVE_MPI
		MPI_Abort(MPI_COMM_WORLD, -2);
#endif 
		return -2;
	}
	
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	return 0;
}




























