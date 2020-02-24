#ifndef __UTILITIES
#define __UTILITIES

#include <vector>
#include <complex>
#include <boost/random.hpp>
#include <valarray>
#include "MPIUtilities.h"

//Input für Matrizen ist Row-Major
//BOOST_STATIC_ASSERT(numeric_limits<double>::is_iec559) for memset ?????

extern "C" {
	double dasum_(int const*, double const*, int const*);
	double dnrm2_(int const*, double const*, int const*);
	double ddot_(const int*, const double*, const int*, const double*, const int*);
	void   dswap_(const int*, double*, int const*, double*, int const*);
	void   daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
	void   dger_(const int*, const int*, double const*, double const*, int const*, double const*, int const*, double*, int const*);
	void   dscal_(int const*, double const*, double*, int const*);
	void   dgemm_(const char*, const char*, int const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
	void   dcopy_(int const*, double const*, int const* , double*, int const*);
	void   dgemv_(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void   dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);
	void   zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
}

namespace Ut {		
	//-----------------------------------------------------------------------------------------------------------------------------------------------------
	
	typedef std::complex<double> complex;
	typedef boost::mt19937 EngineType;
	typedef boost::uniform_real<double> UniformDistribution;
	typedef boost::variate_generator<EngineType&, UniformDistribution> UniformRngType;
	
	//-------------------------------------------------------------------------------------------------------------------------------------------------------
	
	template<class T>
	struct MeasEntry {
		MeasEntry() {};
		MeasEntry(T mean, T error) : mean_(mean), error_(error) {};
		T mean() const { return mean_;};
		T error() const { return error_;};
	private:
		T mean_; T error_;
	};
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------
	
	template<class T>
	void read(std::ifstream& file, T& t) {
		file.read(reinterpret_cast<char*>(&t), sizeof(T));
	};
	
	template<class T>
	void write(std::ofstream& file, T t) {
		file.write(reinterpret_cast<char*>(&t), sizeof(T));
	};

struct Observable {
		Observable() : counter_(0) {};
		
		void add(std::valarray<double> const& val) { 
			if(obs_.size() == 0) 
				obs_.resize(val.size(), .0);  //Is supposed to be initialized to zero anyway, but it's better for my nerves this way
			
			if(obs_.size() != val.size()) 
				throw std::runtime_error("Meas: missmatch in array size !");
			
			++counter_;
			obs_ += val;
		};
		
		void reduce(json_spirit::mValue& jEntry) {
			if(!counter_)
				throw std::runtime_error("Meas: No measurements taken !");
			
			boost::int64_t accCounter;
			std::valarray<double> accObs(mpi::rank() == mpi::master ? obs_.size() : 0);
			
#ifdef HAVE_MPI
			MPI_Reduce(&counter_, mpi::rank() == mpi::master ? &accCounter : 0, 1, MPI_INT64_T, MPI_SUM, mpi::master, MPI_COMM_WORLD);
			MPI_Reduce(&obs_[0], mpi::rank() == mpi::master ? &accObs[0] : 0, obs_.size(), MPI_DOUBLE, MPI_SUM, mpi::master, MPI_COMM_WORLD);			
#else
			accCounter = counter_;
			accObs = obs_;
#endif
			
			if(mpi::rank() == mpi::master) {
				accObs /= accCounter; 
				jEntry = json_spirit::mValue(&accObs[0], &accObs[0] + accObs.size());
			}
		};
	private:
		boost::int64_t counter_; 
		std::valarray<double> obs_;		
	};
	
	void operator<<(Observable& obs, double val) { obs.add(std::valarray<double>(val, 1));};	
	void operator<<(Observable& obs, std::valarray<double> const& val) { obs.add(val);};
	
	typedef std::map<std::string, Observable> Measurements;

	struct Simulation {
		Simulation(std::string inputFolder,std::string outputFolder,std::string name) : name_(name), outputFolder_(outputFolder) {
			{
			    json_spirit::mValue temp;
			    mpi::read_json(inputFolder + name_ + ".json", temp);
			    jParams_ = temp.get_obj();
			}

			jParams_["SEED"] = jParams_.at("SEED").get_int() + mpi::rank();
		};
		
		json_spirit::mObject const& params() { return jParams_;};
		json_spirit::mObject& jobspecs() { return jJobSpecs_;};
		Measurements& meas() { return measurements_;};
		
		void save(boost::int64_t thermalization_sweeps, boost::int64_t measurement_sweeps) {
			boost::int64_t acc_thermalization_sweeps;
			boost::int64_t acc_measurement_sweeps;
			
#ifdef HAVE_MPI
			MPI_Reduce(&thermalization_sweeps, mpi::rank() == mpi::master ? &acc_thermalization_sweeps : 0, 1, MPI_INT64_T, MPI_SUM, mpi::master, MPI_COMM_WORLD);
			MPI_Reduce(&measurement_sweeps, mpi::rank() == mpi::master ? &acc_measurement_sweeps : 0, 1, MPI_INT64_T, MPI_SUM, mpi::master, MPI_COMM_WORLD);
#else
			acc_thermalization_sweeps = thermalization_sweeps;
			acc_measurement_sweeps = measurement_sweeps;
#endif
						
			json_spirit::mObject jMeas;			
			for(std::map<std::string, Observable>::iterator it = measurements_.begin(); it != measurements_.end(); ++it) 
				it->second.reduce(jMeas[it->first]);
				
			if(mpi::rank() == mpi::master) {
				jJobSpecs_["Thermalization Sweeps per Processor"] = acc_thermalization_sweeps/static_cast<double>(mpi::number_of_workers());
				jJobSpecs_["Measurement Sweeps per Processor"] = acc_measurement_sweeps/static_cast<double>(mpi::number_of_workers());
				jJobSpecs_["Number of Processors"] = mpi::number_of_workers();
				
				json_spirit::mObject jParams = jParams_;
				jParams["SEED"] = jParams.at("SEED").get_int() + mpi::number_of_workers(); 
				
			    json_spirit::mObject temp;
				temp["Parameters"] = jParams;
				temp["Job Specifications"] = jJobSpecs_; 
			    temp["Measurements"] = jMeas;

				std::ofstream file((outputFolder_ + name_ + ".meas.json").c_str());				
				json_spirit::write(temp, file, json_spirit::pretty_print | json_spirit::single_line_arrays | json_spirit::remove_trailing_zeros);
				file.close();
			}
#ifdef HAVE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
		std::string outputFolder()
		{
			return outputFolder_;
		}
	private:
		std::string const name_;
		std::string const outputFolder_;
		json_spirit::mObject jParams_;
		json_spirit::mObject jJobSpecs_;
		Measurements measurements_;
	};
	
};
#endif