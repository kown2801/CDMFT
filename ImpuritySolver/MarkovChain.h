#ifndef __MARKOVCHAIN
#define __MARKOVCHAIN

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cassert>
#include <utility>
#include <json_spirit.h>
#include "Utilities.h"
#include "Link.h"
#include "Bath.h"
#include "Trace.h"


namespace Ma {
	struct MarkovChain {
	public:
		static void print_copyright(std::ostream& out) {
			out << "lall, schwaetz, d.h. vorallem dings wenn's z'streng wird." << std::endl << std::endl;
		};
			
		/** 
		* 
		* MarkovChain(json_spirit::mObject const& jNumericalParams, json_spirit::mObject const& jHyb, json_spirit::mArray const& jLink,Ut::Simulation& simulation)
		* 
		* Parameters :	jNumericalParams : storage of all the numerical parameters of the simulation
		*				jHyb : storage of the hybridization informations (read more in README.MD)
		*				jLink : storage of the link informations (read more in README.MD)
		*				simulation : Object used to store the parameters and the measurements throughout the simulation
		* 
		* Prints :	the SEED of the QMC program
		*			any error related to mismatches between the parameters in different files
		*			if the config files do not exist
		*
		* Description: 
		*   The Markovchain constructor is used to initialize the simulation.
		*	If the config file for the current processor is located in the outputFolder, it loads the operators from this file for the starting point
		*	Else the program starts with an empty segment picture
		* 
		*/
		
		MarkovChain(json_spirit::mObject const& jNumericalParams, json_spirit::mObject const& jHyb, json_spirit::mArray const& jLink,Ut::Simulation& simulation) : 
		simulation_(simulation),
		node_(mpi::rank()),
		rng_(jNumericalParams.at("SEED").get_int()),
		urng_(rng_, Ut::UniformDistribution(.0, 1.)),
		beta_(jNumericalParams.at("beta").get_real()),
		probFlip_(jNumericalParams.at("PROBFLIP").get_real()),
		nSite_(jLink.size()),
		link_(jNumericalParams, jHyb, jLink, simulation.meas()),
		trace_(nSite_, static_cast<Tr::Trace*>(0)),
		bath_(2, static_cast<Ba::Bath*>(0)),
		signTrace_(1),
		signBath_(1),
		accSign_(.0),
		probabilityHavingKOperators_(.0, jNumericalParams.at("EOrder").get_int() > 0 ? jNumericalParams.at("EOrder").get_int() : 1),
		acc_(jNumericalParams),
		accChiij_(.0,nSite_*nSite_*acc_.Chi.size()),
		accChiijReal_(0.,accChiij_.size()),
		updateAcc_(2*nSite_, .0),
		updateTot_(2*nSite_, .0),
		updateFlipAcc_(0),
		updateFlipTot_(0) {
			mpi::cout << jNumericalParams.at("SEED").get_int() << std::endl;
			Ut::Measurements& measurements = simulation.meas();
			bath_[0] = new Ba::Bath();
			bath_[1] = new Ba::Bath();
			
			json_spirit::mValue jPreviousConfigValue;
			try{
				mpi::read_json_all_processors(simulation_.outputFolder() + "config_" + boost::lexical_cast<std::string>(node_) + ".json", jPreviousConfigValue);
				json_spirit::mObject const& jPreviousConfig = jPreviousConfigValue.get_obj();
				if(jPreviousConfig.size()) {
					double beta = jPreviousConfig.at("beta").get_real();
					if(beta != beta_) throw std::runtime_error("MarkovChain: missmatch in beta");
					
					int nSite = jPreviousConfig.at("nSite").get_int();
					if(nSite != nSite_) throw std::runtime_error("MarkovChain: missmatch in site number");
					
					for(int i = 0; i < nSite_; ++i) {
						trace_[i] = new Tr::Trace(jNumericalParams, i, measurements, jPreviousConfig);
						
						signTrace_ *= trace_[i]->sign();
						
						for(Tr::Trace::Operators::const_iterator it = trace_[i]->operators(0).begin(); it != trace_[i]->operators(0).end(); ++it) 
							it->type() ? bath_[0]->addDagg(i, it->time(), it->ptr()) : bath_[0]->add(i, it->time(), it->ptr());
						for(Tr::Trace::Operators::const_iterator it = trace_[i]->operators(1).begin(); it != trace_[i]->operators(1).end(); ++it) 
							it->type() ? bath_[1]->addDagg(i, it->time(), it->ptr()) : bath_[1]->add(i, it->time(), it->ptr());
					};
					
					signBath_ *= bath_[0]->rebuild(link_);
					signBath_ *= bath_[1]->rebuild(link_);
					
					int key = jPreviousConfig.at("key").get_int();
					if(key != 5345433) 
						throw std::runtime_error("MarkovChain: error while reading config file.");
				} else 
				{
					for(int i = 0; i < nSite_; ++i) 
					{
						trace_[i] = new Tr::Trace(jNumericalParams, i, measurements, jPreviousConfig);
					}
				}
			}catch(...)
			{
				json_spirit::mObject jPreviousConfig = json_spirit::mObject();
				mpi::cout << "No config file found" << std::endl;
				for(int i = 0; i < nSite_; ++i) 
				{
						trace_[i] = new Tr::Trace(jNumericalParams, i, measurements, json_spirit::mObject());
				}
			}
		}
		
		/** 
		* 
		* void doUpdate()
		*
		* Description: 
		*   Suggests and accepts updates on the segment picture
		*	The possible updates are : 
		*		Spin flipping of a site (completely flips the up and down operators)
		*		Flipping all the operators of two sites
		*		Removing a pair of operators
		*		Inserting a new pair of operators
		* 
		*/
		void doUpdate() { 
			int const site = static_cast<int>(urng_()*nSite_);

			if(urng_() < probFlip_) {
				++updateFlipTot_;
				urng_() < .5 ? flipSpin(site) : flipSite(site, static_cast<int>(urng_()*nSite_));
			} else {				
				int const spin = static_cast<int>(urng_()*2);
				++updateTot_[2*site + spin];

				if(trace_[site]->operators(spin).size())
				    urng_() < .5 ? insert(site, spin) : erase(site, spin);
				else if(urng_() < .5) 
					insert(site, spin);
			};
		};

		/** 
		* 
		* void measure()
		*
		* Description: 
		*	Saves the observables on each sites and in the bath. 
		* 
		*/
		void measure() {
			int const sign = signTrace_*signBath_;
			
			accSign_ += sign;
			
			int k = 0;
			std::vector<std::valarray<Ut::complex> > Chiij_local(nSite_,std::valarray<Ut::complex>(0.,acc_.Chi.size())); //Stores locally the Chiij data coming from the Trace functions
			for(int i = 0; i < nSite_; ++i) {
				Chiij_local[i] = trace_[i]->measure(sign);
				k += (trace_[i]->operators(0).size() + trace_[i]->operators(1).size())/2;
			};
			for(int i = 0;i<nSite_;i++){
				for(int j = 0;j<nSite_;j++){
					accChiij_[std::slice((i*nSite_ + j)*acc_.Chi.size(),acc_.Chi.size(),1)] += Ut::complex(sign,0)*Chiij_local[i]*Chiij_local[j].apply(std::conj);
				}
			}
			
			if(k < static_cast<int>(probabilityHavingKOperators_.size())) probabilityHavingKOperators_[k] += sign;
			
			link_.measure(sign, bath_[0]->begin(), bath_[0]->end());
			link_.measure(sign, bath_[1]->begin(), bath_[1]->end());
		};

		/** 
		* 
		* void store(Ut::Measurements& measurements, int measurementsFromLastStore)
		* 
		* Parameters :	measurements : variable used to store the measurements, used for output
		*				measurementsFromLastStore : Number of measurements done simce the last time we stored some measurements
		*
		* Description: 
		*   Stores the measurements in the handler.
		*	Resets the measurements for the next round of measurements
		* 
		*/
		void store(Ut::Measurements& measurements, int measurementsFromLastStore) {
			measurements["Sign"] << accSign_/measurementsFromLastStore; accSign_ = 0;
			
			probabilityHavingKOperators_ /= measurementsFromLastStore;
		    measurements["pK"] << probabilityHavingKOperators_;
			probabilityHavingKOperators_ = .0;

			//The processing of the spin susceptibility is a little bit more complicated
			std::valarray<Ut::complex> temp = accChiij_[std::slice(0,nSite_*nSite_,acc_.Chi.size())];
			temp /= beta_;
			accChiij_[std::slice(0,nSite_*nSite_,acc_.Chi.size())] = temp;
			//We need to get an array of doubles
			for(unsigned int n=0;n<accChiij_.size();n++){
				accChiijReal_[n] = accChiij_[n].real();
			}
			//As stated in Trace.h, the result is not exactly good for the same site susceptibility when there is no particles of a certain spin on the sites. We correct that here
			for(int n=0;n<nSite_;n++){
				accChiijReal_[std::slice((n*nSite_ + n)*acc_.Chi.size(),acc_.Chi.size(),1)] = trace_[n]->getChi();
			}
			for(unsigned int n=1;n<acc_.Chi.size();n++){
				std::valarray<double> temp = accChiijReal_[std::slice(n,nSite_*nSite_,acc_.Chi.size())];
			        temp *= beta_/((2*n*M_PI)*(2*n*M_PI));
				accChiijReal_[std::slice(n,nSite_*nSite_,acc_.Chi.size())] = temp;
			}
			accChiijReal_ /= measurementsFromLastStore;
			measurements["Chiij"] << accChiijReal_;
			accChiij_ = .0;
			//End of processing of the spin susceptibility
			
			for(int i = 0; i < nSite_; ++i) trace_[i]->store(measurements, i, acc_, measurementsFromLastStore);
			
			acc_.N /= nSite_;
			acc_.Sz /= nSite_;
			acc_.D /= nSite_;
			acc_.Chi /= nSite_;
			
			measurements["k"] << acc_.k;
			measurements["N"] << acc_.N;
			measurements["D"] << acc_.D;
			measurements["Sz"] << acc_.Sz;
			measurements["Chi0"] << acc_.Chi[0];
			
			if(acc_.Chi.size() > 1) 
				measurements["Chi"] << acc_.Chi;
			
			acc_.k = .0;
			acc_.N = .0;
			acc_.Sz = .0;
			acc_.D = .0;
			acc_.Chi = .0;

			link_.store(measurements, measurementsFromLastStore);
		};
		/** 
		* 
		* void cleanUpdate()
		*
		* Description: 
		*	Due to numerical error, we need to rebuild the entire bath matrix from time to time
		*	Because this is an expensive operation, we don't do this operation too often
		* 
		*/
		void cleanUpdate() { 
			signBath_ = 1;
			signBath_ *= bath_[0]->rebuild(link_);
			signBath_ *= bath_[1]->rebuild(link_);
		};
		/** 
		* 
		* ~MarkovChain()
		* 
		* Description: 
		*	Saves the current configuration in the config files to be able to resume the simulation. 
		* 
		*/
		~MarkovChain() {
			for(int site = 0; site < nSite_; ++site) {
				mpi::cout << site << "Up:  " << updateAcc_[2*site]/static_cast<double>(updateTot_[2*site]) << std::endl;
				mpi::cout << site << "Down:  " << updateAcc_[2*site + 1]/static_cast<double>(updateTot_[2*site + 1]) << std::endl;
			};
			
			mpi::cout << "Flip:  " << updateFlipAcc_/static_cast<double>(updateFlipTot_) << std::endl;
			
			delete bath_[1];
			delete bath_[0];	
			
			json_spirit::mObject jConfig;
			int key = 5345433; 
			jConfig["beta"] = beta_;
			jConfig["nSite"] = nSite_;
			jConfig["key"] = key;
			for(int site = 0; site < nSite_; ++site) {
				trace_[site]->saveConfig(jConfig,site);
				delete trace_[site];
			};

			std::ofstream file((simulation_.outputFolder() + "config_" + boost::lexical_cast<std::string>(mpi::rank()) + ".json").c_str());			
			json_spirit::write(jConfig, file, json_spirit::pretty_print | json_spirit::single_line_arrays | json_spirit::remove_trailing_zeros);			
			file.close();
		};
	private:
		Ut::Simulation& simulation_;
		int const node_;
		Ut::EngineType rng_;
		Ut::UniformRngType urng_;
		
		double const beta_;
		double const probFlip_;
		int const nSite_;
		
		Link::Link link_;
		std::vector<Tr::Trace*> trace_;
		std::vector<Ba::Bath*> bath_;
		
		int signTrace_;
		int signBath_;
		
		double accSign_;
		std::valarray<double> probabilityHavingKOperators_;
		Tr::Meas acc_;
		std::valarray<Ut::complex> accChiij_;//This valarray has a weird structure. It contains Chiij for all bosonic Matsubara frequencies. Sorted first by sites indices and then Matsubara frequencies.
											//accChiij_[(i*nSite_ + j)*acc_.chi.size() + n] give Chi(i,j,i\oemga_n)
		std::valarray<double> accChiijReal_;//As we can only save real numbers, we need a second array to store it
		
		std::vector<double> updateAcc_;
		std::vector<double> updateTot_;
		int updateFlipAcc_;
		int updateFlipTot_;
			
		/** 
		* 
		* void insert(int site, int spin)
		* 
		* Parameters :	site : site number
		*				spin : spin
		*
		* Description:
		*   Tries to insert a new pair of operators in the imaginary time line
		*	Accepts or reject the insertion according to the Metropolis Hastings acceptation rate.
		* 
		*/
		void insert(int site, int spin) {
			Tr::Trace& trace = *trace_[site];
			Ba::Bath& bath = *bath_[spin];
				
			double prob = .0;
			prob += trace.insert(spin, urng_);
			prob += bath.insert(site, trace.op().time(), trace.op().ptr(), trace.opDagg().time(), trace.opDagg().ptr(), link_);

			if(std::log(1. - urng_()) < prob) { 
				signTrace_ *= trace_[site]->acceptInsert(spin);
			    signBath_ *= bath.acceptInsert();

				++updateAcc_[2*site + spin];
				
				return;
			}
			
			trace.rejectInsert(spin);
		};
		/** 
		* 
		* void erase(int site, int spin)
		* 
		* Parameters :	site : site number
		*				spin : spin
		*
		* Description: 
		*   Tries to erase a a pair of operators in the imaginary time line.
		*	Accepts or reject the erase according to the Metropolis Hastings acceptation rate.
		* 
		*/
		void erase(int site, int spin) {
			Tr::Trace& trace = *trace_[site];
			Ba::Bath& bath = *bath_[spin];

			double prob = .0;
			prob += trace.erase(spin, urng_);			
			prob += bath.erase(site, trace.op().time(), trace.opDagg().time());
			
			if(std::log(1. - urng_()) < prob) { 
				signTrace_ *= trace_[site]->acceptErase(spin);
				signBath_ *= bath.acceptErase();

				++updateAcc_[2*site + spin];
							
				return;
			}
			
			trace.rejectErase(spin);
		};
		
		/** 
		* 
		* void flipSpin(int site)
		* 
		* Parameters :	site : site number
		*
		* Description: 
		*  	Tries to flip the spins of one site.
		*	Accepts or reject the flip according to the result of the tryFlip() function.
		* 
		*/
		void flipSpin(int site) {
		    trace_[site]->flip();
			//If the flip is not accepted, we swap again
		    if(!tryFlip()) trace_[site]->flip();
		};
		
		/** 
		* 
		* void flipSite(int siteA, int siteB)
		* 
		* Parameters :	siteA : first site number
		*				siteB : second site number
		*
		* Description: 
		*  	Tries to flip all the operators of two sites.
		*	Accepts or reject the flip according to the result of the tryFlip() function.
		* 
		*/
		void flipSite(int siteA, int siteB) {
			std::swap(trace_[siteA], trace_[siteB]);
			//If the flip is not accepted, we swap again
			if(!tryFlip()) std::swap(trace_[siteA], trace_[siteB]);
		};
		
		/** 
		* 
		* int tryFlip()
		* 
		* Return Value : returns 1 if the flip is accepted
		*				 returns 0 if the flip is not accepted  
		* Description: 
		*	Fills new bath with the flipped operators
		*  	Computes the Metropolis Hastings Acceptation rate for the flip using the old and the new baths.
		*	In case the flip is accepted : 
		*		the new bath matrices are updated with their new values (rebuild)
		*		the signs are correctly set
		*		the simulation can continue without anything else
		*/
		int tryFlip() {
			std::vector<Ba::Bath*> newBath(2);
			newBath[0] = new Ba::Bath();
			newBath[1] = new Ba::Bath();
			
			for(int i = 0; i < nSite_; ++i) {
				Tr::Trace& trace = *trace_[i];
				
				for(Tr::Trace::Operators::const_iterator it = trace.operators(0).begin(); it != trace.operators(0).end(); ++it) 
					it->type() ? newBath[0]->addDagg(i, it->time(), it->ptr()) : newBath[0]->add(i, it->time(), it->ptr());
				
				for(Tr::Trace::Operators::const_iterator it = trace.operators(1).begin(); it != trace.operators(1).end(); ++it) 
					it->type() ? newBath[1]->addDagg(i, it->time(), it->ptr()) : newBath[1]->add(i, it->time(), it->ptr());
			}
		
			newBath[0]->rebuild(link_);
			newBath[1]->rebuild(link_);
		
			if(urng_() < std::abs(newBath[0]->det()/bath_[0]->det()*newBath[1]->det()/bath_[1]->det())) {
				//We accept the new bath and get rid of the old bath
				signTrace_ = 1;
				for(std::vector<Tr::Trace*>::iterator it = trace_.begin(); it != trace_.end(); ++it)
					signTrace_ *= (*it)->sign();
				
				delete bath_[0]; bath_[0] = newBath[0];
				delete bath_[1]; bath_[1] = newBath[1];
				
				signBath_ = 1;
				signBath_ *= bath_[0]->det() > .0 ? 1 : -1;
				signBath_ *= bath_[1]->det() > .0 ? 1 : -1;
				
				++updateFlipAcc_;
				
				return 1;
			};

			//the new bath is not accepted
			delete newBath[0];
			delete newBath[1];	
			
			return 0;
		};
	};
};

#endif




















