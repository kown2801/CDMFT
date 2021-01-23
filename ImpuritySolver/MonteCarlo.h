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
	
	void MonteCarlo(Ma::MarkovChain& markovChain, Ut::Simulation& simulation) //Lauching the program
	{	
		std::time_t time;
		Timer timer;
		
		int64_t thermalization_sweeps = 0;
		int64_t measurement_sweeps = 0;
		int64_t measurementsFromLastSample = 0;

		json_spirit::mObject const& jParams = simulation.params();		
		int64_t const clean_every_sweep = jParams.at("CLEAN_EVERY_SWEEP").get_int64();
		int64_t const sample_every_sweep = jParams.at("SAMPLE_EVERY_SWEEP").get_int64();
		int64_t const store_every_sample = jParams.at("STORE_EVERY_SAMPLE").get_int64();
		mpi::cout << "Start thermalization at " << std::ctime(&(time = std::time(NULL))) << std::flush; 
		mpi::cout = mpi::one;
		mpi::cout << "We go for " << 60.*jParams.at("THERMALIZATION_TIME").get_real() << " seconds of MonteCarlo thermalization" << std::endl;
		mpi::cout  = mpi::every;	
		timer.start(60.*jParams.at("THERMALIZATION_TIME").get_real());
		for(; 1; ) { 
			++thermalization_sweeps;
			markovChain.doUpdate();
			
			if(thermalization_sweeps % clean_every_sweep == 0) {
				markovChain.cleanUpdate();
			}
			
			if(thermalization_sweeps % sample_every_sweep == 0)
			{
				if(timer.end()) break;
			}
		}
		mpi::cout << "Start measurements at " << std::ctime(&(time = std::time(NULL))) << std::flush;
		mpi::cout = mpi::one;
		mpi::cout << "We go for " << 60.*jParams.at("MEASUREMENT_TIME").get_real() << " seconds of MonteCarlo simulation" << std::endl;
		mpi::cout = mpi::every;
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
					markovChain.measure(simulation.meas(), store_every_sample);
					if(timer.end()) break;					
				}
			}
		}
		mpi::cout << "Start saving simulation at " << std::ctime(&(time = std::time(NULL))) << std::flush;
		simulation.save(thermalization_sweeps, measurement_sweeps);
	};
};

#endif











