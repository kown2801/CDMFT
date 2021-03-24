#ifndef __UTILITIES
#define __UTILITIES

#include <vector>
#include <complex>
#include <random>
#include <valarray>
#include "MPIUtilities.h"
#include <iterator>
#include "nlohmann_json.hpp"
using json=nlohmann::json;

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
	typedef std::mt19937 EngineType;
	typedef std::uniform_real_distribution<double> UniformDistribution;
	typedef std::function<double()> UniformRngType;
	
	//-------------------------------------------------------------------------------------------------------------------------------------------------------
	
	template<class T>
	struct MeasEntry {
		MeasEntry() {};
		MeasEntry(T mean, T error) : mean_(mean), error_(error) {};
		T binning_mean() const { return mean_;};
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
		/** 
		* 
		* void add(std::valarray<double> const& val)
		* 
		* Parameters :	val : input array to store in the measurements
		* 
		* Description: 
		*   Stores the measurements. It is a bit complicated because we want to do a binning analysis in order to get an accurate approximation of the error
		* 	We store the sum of measurements but also the sum of the squares of the measurements.
		*	This is a direct copy of the ALPS library
		*/
		void add(std::valarray<double> const& x){
  		    // set sizes if starting additions
  		    if(counter_==0)
  		    {
  		      last_bin_.resize(1);
  		      sum_.resize(1);
  		      sum2_.resize(1);
  		      bin_entries_.resize(1);
  		      last_bin_[0].resize(x.size());
  		      sum_[0].resize(x.size());
  		      sum2_[0].resize(x.size());
  		    }
		  	  		   // store x, x^2
  		    last_bin_[0]=x;
  		    sum_[0]+=x;
  		    sum2_[0]+= x*x;

		    uint64_t i=counter_;
  		    counter_++;
  		    bin_entries_[0]++;
  		    uint64_t binlen=1;
  		    std::size_t bin=0;
		  	  		   // binning
  		    do{
  		        if(i&1){
  		            // a bin is filled
	  		            binlen*=2;
	  		            bin++;
	  		            if(bin>=last_bin_.size()){
	  		                last_bin_.resize((std::max)(bin+1,last_bin_.size()));
	  		                sum_.resize((std::max)(bin+1, sum_.size()));
	  		                sum2_.resize((std::max)(bin+1,sum2_.size()));
	  		                bin_entries_.resize((std::max)(bin+1,bin_entries_.size()));
  
  			  	  		    last_bin_[bin].resize(x.size());
  	  		                sum_[bin].resize(x.size());
  	  		                sum2_[bin].resize(x.size());
	  		            }

			  	  		std::valarray<double>  x1=(sum_[0]-sum_[bin]);
	  		            x1/=double(binlen);
			  	  		std::valarray<double>  y1 = x1*x1;
			  	  		last_bin_[bin]=x1;
	  		            sum2_[bin] += y1;
	  		            sum_[bin] = sum_[0];
	  		            bin_entries_[bin]++;
  		          }
  		        else{
  		          break;
  		        }
  		    } while (i>>=1);
		}
		/* START of helper functions to conmpute the error correctly */
		uint64_t binning_count() const {return counter_;} // number of measurements performed
		std::valarray<double> binsum(std::size_t i) const
		{
		  return sum_[i]/double(1ll<<i);
		}
		std::valarray<double> binsquared_sum(std::size_t i) const
		{
		  std::valarray<double> retval(sum2_[i]);
		  return retval;
		}
		uint64_t binsize(std::size_t i) const //Returns the number of entry in the bin
		{
			return bin_entries_[i];
		}
		uint32_t binning_depth() const //Returns the maximum binning level we study.
		{
		  return ( int(sum_.size())< 1 ) ? 1 : int(sum_.size());
		}
		std::valarray<double> binning_sum() const
		{
		  return sum_[0];
		}
		/* End of helper functions */

		/** 
		* 
		* void reduce(json& jMeas,json& jErrors)
		* 
		* Parameters :	jEntry : storage point of the computed observables for printing in the ouput file
		*				jMeas : storage point of the errors of the observables for printing in the ouput file
		* 
		* Description: 
		*   Takes all the measured quantities from all running simulations in the MPIWorld together to a single value for the output file.
		* 	To do this, we sum all the measurements for this observable over the MPIWorld and divide by the total number of measurements
		*	We then add the entry to the output file. This is fairly easy to do here.
		* 	Then we want to compute the error. This is a bit more difficult.
		* 	We do a binning analysis using the conserved quantities computed in add on the single processors and also between them.
		*	We then save the error in the jErrors object. 
		* 	If you want to know more, the process is detailed in the function. 
		*/
		void reduce(json& jMean,json& jErrors) {
			
			if(counter_ < 2)
				throw std::runtime_error("Meas: Not enough measurements taken !");
			//First we need to accumulate the counter in order to know how much measurements we have in total
			uint64_t accCounter = 0;
			mpi::reduce(counter_,accCounter);

			//Now we accumulate the observable sums
			std::valarray<double> binningSum = binning_sum();
			std::valarray<double> accBinningMean(mpi::rank() == mpi::master ? binningSum.size() : 0);
			mpi::reduce(binningSum,accBinningMean);
			//Then we compute and save the mean only on the first processor
			if(mpi::rank() == mpi::master) {
				accBinningMean/=accCounter;
				jMean = accBinningMean;
			}
			
			//Now we get the error. This operation is a bit trickier
			//First we get the error from the simulations when considered independently (what Alps does)
				//First we need to get the maximum binning_depth across processors
			uint64_t binningDepth = binning_depth();
			uint64_t maxBinningDepth;
			mpi::getMax(binningDepth,maxBinningDepth);

			std::vector<double> accBinningErrs;
				//Then we compute the error for all binning depths.
				//We do this by first merging the data from all processors at each binning level (add the sums and the counts)
				//And then computing the error using the un-biased variance 
			for(uint32_t i=0;i< maxBinningDepth;i++)
			{
				uint64_t nbTerms(0);
				uint64_t accNbTerms(0);
				//If the current binning level is above the maximum binning level for this processors, we do nothing
				if(binningDepth > i){
					nbTerms = binsize(i);
				}

				//Getting the total number of terms for this bin number
				mpi::reduce(nbTerms,accNbTerms);

				//Processing the data (adding the plain sums and the squared sums across processors)
				std::valarray<double> binningSquaredSum(binning_sum().size());
				std::valarray<double> binningSum(binning_sum().size());
				std::valarray<double> accBinningSquaredSum(mpi::rank() == mpi::master ? binningSquaredSum.size() : 0);
				std::valarray<double> accBinningSum(mpi::rank() == mpi::master ? binningSum.size() : 0);
				if(binningDepth > i){
					binningSquaredSum = binsquared_sum(i);
					binningSum = binsum(i);
				}
				mpi::reduce(binningSquaredSum,accBinningSquaredSum);//getting the total squared sum for the binning level
				mpi::reduce(binningSum,accBinningSum);//getting the total sum for the binning level

				//We compute the error on the first processors and save it to the accBinningErrs array
				if(mpi::rank() == mpi::master) {
					std::valarray<double> accBinningErr = std::sqrt((accBinningSquaredSum/double(accNbTerms-1) - accBinningSum*accBinningSum/(double(accNbTerms-1)*double(accNbTerms)))/double(accNbTerms));
					accBinningErrs.insert(accBinningErrs.end(),std::begin(accBinningErr),std::end(accBinningErr));
				}
			}

			//Then we want to get the error (and maybe some error convergence) thanks to our multiple processors
				//We need a list of all means and measurement counts
				//We also need to gather all the data in a one dimensional vector and then put in a vector of valarray (needed for MPI_Gather)
			std::valarray<double> allBinningSumsGather(static_cast<double>(0),mpi::number_of_workers()*binningSum.size());
#ifdef HAVE_MPI
			MPI_Gather(&binningSum[0],binningSum.size(),MPI_DOUBLE,&allBinningSumsGather[0],binningSum.size(), MPI_DOUBLE,0,MPI_COMM_WORLD);	
#else
			allBinningSumsGather = binningSum;
#endif
			std::vector<uint64_t> allBinningCounts(mpi::number_of_workers(),0);
#ifdef HAVE_MPI
			MPI_Gather(&counter_,1,MPI_DOUBLE,&allBinningCounts[0],1, MPI_DOUBLE,0,MPI_COMM_WORLD);
#else
			allBinningCounts[0] = counter_;
#endif
				//We do a binning analysis on the collected data (is it really necessary when the processors are supposed to have independent simulations ?)
			if(mpi::rank() == mpi::master) {
				//We need to put the Sums into a vector of valarray : 	
				std::vector<std::valarray<double> > allBinningSums;
				for(int i=0;i<mpi::number_of_workers();i++){
					allBinningSums.push_back(allBinningSumsGather[std::slice(i*binningSum.size(),binningSum.size(),1)]);
				}
				//We compute the statistics on the present arrays
				while(allBinningSums.size() > 1){
					uint64_t N_total = std::accumulate(allBinningCounts.begin(), allBinningCounts.end() , 0);
					uint64_t N_stage = allBinningCounts.size();
					std::valarray<double> current_total_sum(static_cast<double>(0),binningSum.size());
					current_total_sum = std::accumulate(allBinningSums.begin(), allBinningSums.end() , current_total_sum);
					std::valarray<double> error(static_cast<double>(0),binningSum.size());
					for(uint64_t i=0;i<N_stage;i++){
						error += allBinningSums[i]*allBinningSums[i]/allBinningCounts[i]/N_total;
					}
					error -= (current_total_sum/(double)N_total)*(current_total_sum/(double)N_total);
					error = std::sqrt((double)(N_stage/(N_stage - 1))*error/N_stage);
					//We save the current data status
					accBinningErrs.insert(accBinningErrs.end(),std::begin(error),std::end(error));

					//We reduce the data size by 2 by regrouping data sets 2 by 2
					for(uint64_t i=0;2*i+1<N_stage;i++){
						allBinningCounts[i] = allBinningCounts[2*i] + allBinningCounts[2*i+1];
						allBinningSums[i] = allBinningSums[2*i] + allBinningSums[2*i+1];
					}
					allBinningSums.resize(N_stage/2);
					allBinningCounts.resize(N_stage/2);
				}
			}
			if(mpi::rank() == mpi::master) {
				/* In order to save the whole error graphs if you want to do some error analysis, use the following line */
					//jErrors = accBinningErrs;
				/* Else we need for each vector component, only the maximum error */
					///*
					std::valarray<double> binningErrs(accBinningMean.size());
					//We need to have a valarray to be able to perform the max operation on a slice
					std::valarray<double> accBinningErrsValArray(accBinningErrs.data(), accBinningErrs.size());
					uint32_t n_errs = accBinningErrs.size()/binningErrs.size();
					for(uint32_t i = 0 ; i<binningErrs.size();i++){
						binningErrs[i] = static_cast<std::valarray<double> >(accBinningErrsValArray[std::slice(i,n_errs,binningErrs.size())]).max();
					}
					jErrors = binningErrs;
					//*/
			}


		};
	private:
		uint64_t counter_; 
		//Binning implementation
		std::vector<std::valarray< double > > sum_; // sum of measurements in the bin
		std::vector<std::valarray< double > > sum2_; // sum of the squares
		std::vector<uint64_t> bin_entries_; // number of measurements
		std::vector<std::valarray< double > > last_bin_; // the last value measured
	};
	
	void operator<<(Observable& obs, double val) { obs.add(std::valarray<double>(val, 1));};	
	void operator<<(Observable& obs, std::valarray<double> const& val) { obs.add(val);};
	
	typedef std::map<std::string, Observable> Measurements;
	/** 
	* 
	* struct Simulation
	* 
	* Description: 
	*   Object for handling the simulation.
	* 	It allows to : 
	*		Store the input parameters
	* 		Output all the parameters and measurements to the output file
	* 
	*/
	struct Simulation {
		Simulation(std::string inputFolder,std::string outputFolder,std::string name) : name_(name), outputFolder_(outputFolder) {
			{
			    mpi::read_json(inputFolder + name_ + ".json", jParams_);
			}

			jParams_["SEED"] = jParams_["SEED"].get<double>() + mpi::rank();
		};
		
		json const& params() { return jParams_;};
		json& jobspecs() { return jJobSpecs_;};
		Measurements& meas() { return measurements_;};
		/** 
		* 
		* void save(uint64_t thermalization_sweeps, uint64_t measurement_sweeps)
		* 
		* Parameters :	thermalization_sweeps : number of thermalization sweeps
						measurement_sweeps : number of thermalization sweeps (just kidding, number of measurement sweeps)
		* 
		* Description: 
		*   Saves all the measurements and parameters to file.
		* 	First rounds up all the data from all processors.
		*	Then, writing the file is done only on processor 1 to avoid collisions and nonesensical data.
		* 
		*/
		void save(uint64_t thermalization_sweeps,uint64_t measurement_sweeps) {
			int64_t acc_thermalization_sweeps;
			uint64_t acc_measurement_sweeps;
			
#ifdef HAVE_MPI
			MPI_Reduce(&thermalization_sweeps, mpi::rank() == mpi::master ? &acc_thermalization_sweeps : 0, 1, MPI_INT64_T, MPI_SUM, mpi::master, MPI_COMM_WORLD);
			MPI_Reduce(&measurement_sweeps, mpi::rank() == mpi::master ? &acc_measurement_sweeps : 0, 1, MPI_INT64_T, MPI_SUM, mpi::master, MPI_COMM_WORLD);
#else
			acc_thermalization_sweeps = thermalization_sweeps;
			acc_measurement_sweeps = measurement_sweeps;
#endif
						
			json jMeas;		
			json jErrors;
			for(std::map<std::string, Observable>::iterator it = measurements_.begin(); it != measurements_.end(); ++it) 
				it->second.reduce(jMeas[it->first],jErrors[it->first]);

				
			if(mpi::rank() == mpi::master) {
				jJobSpecs_["Thermalization Sweeps per Processor"] = acc_thermalization_sweeps/static_cast<double>(mpi::number_of_workers());
				jJobSpecs_["Measurement Sweeps per Processor"] = acc_measurement_sweeps/static_cast<double>(mpi::number_of_workers());
				jJobSpecs_["Number of Processors"] = mpi::number_of_workers();
				
				json jParams = jParams_;
				jParams["SEED"] = jParams["SEED"].get<double>() + mpi::number_of_workers(); 
				
			    json temp;
				temp["Parameters"] = jParams;
				temp["Job Specifications"] = jJobSpecs_; 
			    temp["Measurements"] = jMeas;
			    temp["Errors"] = jErrors;
		
				IO::writeJsonToFile(outputFolder_ + name_ + ".meas.json", temp);
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
		json jParams_;
		json jJobSpecs_;
		Measurements measurements_;
	};
	
};
#endif
