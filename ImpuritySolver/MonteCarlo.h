#ifndef __MONTECARLO
#define __MONTECARLO

#include <ctime>
#include "MarkovChain.h"

#ifdef HAVE_CHRONO
#include <ratio>
#include <chrono>
#endif

namespace MC {
	
#ifdef HAVE_CHRONO
	
	struct Timer {
		void start(double duration) {
			duration_ = duration;
			start_ = std::chrono::steady_clock::now();
		};
		bool end() {
			return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() > duration_;
		};
	private:
		double duration_;
		std::chrono::steady_clock::time_point start_;
	};
	
#else
	
	struct Timer {
		void start(double duration) {
			duration_ = duration;
			start_ = std::time(NULL);
		};
		bool end() {
			return std::difftime(std::time(NULL), start_) > duration_;
		};
	private:
		double duration_;
		std::time_t start_;
	};
	
#endif
	


	/** 
	* 
	* void MonteCarlo(Ma::MarkovChain& markovchain, Ut::Simulation& simulation)
	* 
	* Parameters :	markovchain : initialized markovchain, the one that does all the job
	*				simulation : Object used to store the parameters and the measurements throughout the simulation
	* 
	* Return Value : Nothing 
	*
	* Prints : the times at which the thermalization starts, the measurements starts, the simulation saves, the simulation finishes
	* 
	* Description: 
	* 
	*    The MonteCarlo Function is the engine of the QMC Simulation
	*    It allows to monitor the number of sweeps, measurements and storing points throughout the simulation.
	* 
	*/
	void MonteCarlo(Ma::MarkovChain& markovChain, Ut::Simulation& simulation) //Lauching the program
	{	
		std::time_t time;
		Timer timer;
		
		boost::int64_t thermalization_sweeps = 0;
		boost::int64_t measurement_sweeps = 0;
		boost::int64_t measurementsFromLastSample = 0;

		json_spirit::mObject const& jParams = simulation.params();		
		boost::int64_t const clean_every_sweep = jParams.at("CLEAN_EVERY_SWEEP").get_int64();
		boost::int64_t const sample_every_sweep = jParams.at("SAMPLE_EVERY_SWEEP").get_int64();
		boost::int64_t const store_every_sample = jParams.at("STORE_EVERY_SAMPLE").get_int64();
		mpi::cout << "Start thermalization at " << std::ctime(&(time = std::time(NULL))) << std::flush; 
		
		timer.start(60.*jParams.at("THERMALIZATION_TIME").get_real());
		for(; 1; ) { 
			++thermalization_sweeps;
			
			markovChain.doUpdate();
			
			if(thermalization_sweeps % clean_every_sweep == 0) 
				markovChain.cleanUpdate();
			
			if(thermalization_sweeps % sample_every_sweep == 0)
				if(timer.end()) break;
		}
		
		mpi::cout << "Start measurements at " << std::ctime(&(time = std::time(NULL))) << std::flush;
		
		timer.start(60.*jParams.at("MEASUREMENT_TIME").get_real());
		for(; 1; ) {
			++measurement_sweeps;
			
			markovChain.doUpdate();
			
			if(measurement_sweeps % clean_every_sweep == 0) 
				markovChain.cleanUpdate();
			
			if(measurement_sweeps % sample_every_sweep == 0) {
				markovChain.measure();
				measurementsFromLastSample++;
				
				if(measurementsFromLastSample % store_every_sample == 0) {
					markovChain.store(simulation.meas(), store_every_sample);	
					measurementsFromLastSample = 0;
					if(timer.end()) break;					
				}
			}
		}
		/*
		Saving the data about hte job, already implemented inside the code, maybe change that
		markovChain.jobspecs(simulation.jobspecs());
		*/


		mpi::cout = mpi::every;
		mpi::cout << "Start saving simulation at " << std::ctime(&(time = std::time(NULL))) << std::flush;
		
		simulation.save(thermalization_sweeps, measurement_sweeps);
	};
};

#endif











