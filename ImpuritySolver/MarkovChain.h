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
#include <nlohmann/json.hpp>
#include "Utilities.h"
#include "Link.h"
#include "Bath.h"
#include "Trace.h"

//check detailed balance for k=0 

namespace Ma {
	struct MarkovChain {
	public:
		static void print_copyright(std::ostream& out) {
			out << "lall, schwaetz, d.h. vorallem dings wenn's z'streng wird." << std::endl << std::endl;
		};
		
		MarkovChain(json const& jNumericalParams, json  const& jHyb, json const& jLink, Ut::Simulation& simulation) :
		simulation_(simulation),
		node_(mpi::rank()),
		rng_(jNumericalParams["SEED"].get<double>()),
		ud_(0.,1.),
		urng_([this](){return ud_(rng_);}),
		beta_(jNumericalParams["beta"]),
		probFlip_(jNumericalParams["PROBFLIP"]),
		nSite_(jLink.size()/2),
		link_(jNumericalParams, jHyb, jLink, simulation.meas()),
		trace_(nSite_, static_cast<Tr::Trace*>(0)),
		bath_(new Ba::Bath()),
		signTrace_(1),
		signBath_(1),
		accSign_(.0),
		pK_(.0, jNumericalParams["EOrder"]),
		acc_(jNumericalParams),
		accChiij_(.0,nSite_*nSite_*acc_.Chi.size()),
		accChiijReal_(0.,accChiij_.size()),
		updateAcc_(2*nSite_, .0),
		updateTot_(2*nSite_, .0),
		updateFlipAcc_(0),
		updateFlipTot_(0) {
			Ut::Measurements& measurements = simulation.meas();
			std::cout << jNumericalParams["SEED"] << std::endl;
            
			json jPreviousConfig;


			try{
				mpi::read_json_all_processors(simulation_.outputFolder() + "config_" + std::to_string(node_) + ".json", jPreviousConfig);
				if(jPreviousConfig.size()) {
					double beta = jPreviousConfig["beta"];
					if(beta != beta_) throw std::runtime_error("MarkovChain: missmatch in beta");
					
					int nSite = jPreviousConfig["nSite"];
					if(nSite != nSite_) throw std::runtime_error("MarkovChain: missmatch in site number");
					
					for(int site = 0; site < nSite_; ++site) {
						trace_[site] = new Tr::Trace(jNumericalParams, site, measurements, jPreviousConfig);
						
						signTrace_ *= trace_[site]->sign();
						
						for(int spin = 0; spin < 2; ++spin)
						    for(Tr::Trace::Operators::const_iterator it = trace_[site]->operators(spin).begin(); it != trace_[site]->operators(spin).end(); ++it) 
							    it->type() ? bath_->addDagg(site, spin, it->time(), it->ptr()) : bath_->add(site, spin, it->time(), it->ptr());
					};
					
					signBath_ *= bath_->rebuild(link_);
					
					int key = jPreviousConfig["key"];
					if(key != 5345433) 
						throw std::runtime_error("MarkovChain: error while reading config file.");
				} else 
				{
					for(int site = 0; site < nSite_; ++site) 
					{
						trace_[site] = new Tr::Trace(jNumericalParams, site, measurements, jPreviousConfig);
					}
				}
			}catch(...)
			{
				json jPreviousConfig;
				mpi::cout << "No config file found" << std::endl;
				for(int site = 0; site < nSite_; ++site) 
				{
						trace_[site] = new Tr::Trace(jNumericalParams, site, measurements, jPreviousConfig);
				}
			}
		}
		
		void doUpdate() { 
			int const site = static_cast<int>(urng_()*nSite_);

			if(urng_() < probFlip_) {
				++updateFlipTot_;
				urng_() < .5 ? flipSpin(site) : flipSite(site, static_cast<int>(urng_()*(nSite_-1)));
			} else {				
				int const spin = static_cast<int>(urng_()*2);
				++updateTot_[2*site + spin];

				if(trace_[site]->operators(spin).size())
				    urng_() < .5 ? insert(site, spin) : erase(site, spin);
				else if(urng_() < .5) 
					insert(site, spin);
			};
		};
		
		void measure() {
			int const sign = signTrace_*signBath_;
			
			accSign_ += sign;
			
			//std::cout << sign << std::endl;
			
			int k = 0;
			std::vector<std::valarray<Ut::complex> > Chiij_local(nSite_,std::valarray<Ut::complex>(0.,acc_.Chi.size()));
			for(int site = 0; site < nSite_; ++site) {
				Chiij_local[site] = trace_[site]->measure(sign);
				k += (trace_[site]->operators(0).size() + trace_[site]->operators(1).size())/2;
			};
			for(int i = 0;i<nSite_;i++){
				for(int j = 0;j<nSite_;j++){
					accChiij_[std::slice((i*nSite_ + j)*acc_.Chi.size(),acc_.Chi.size(),1)] += Ut::complex(sign,0)*Chiij_local[i]*Chiij_local[j].apply(std::conj);
				}
			}
			
			pK_[k] += sign;
			//std::cout << k << std::endl;
			
			link_.measure(sign, bath_->begin(), bath_->end());
		};
		
		void measure(Ut::Measurements& measurements, int NAlpsMeas) {
			//std::cout << accSign_/NAlpsMeas << std::endl;
			
			measurements["Sign"] << accSign_/NAlpsMeas; accSign_ = 0;
			
			pK_ /= NAlpsMeas;
		    	measurements["pK"] << pK_;
			pK_ = .0;

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
			accChiijReal_ /= NAlpsMeas;
			measurements["Chiij"] << accChiijReal_;
			accChiij_ = .0;
			//End of processing of the spin susceptibility

			for(int site = 0; site < nSite_; ++site) trace_[site]->measure(measurements, site, acc_, NAlpsMeas);
			
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

			link_.measure(measurements, NAlpsMeas);	
		};
		
		void cleanUpdate() { 
			signBath_ = bath_->rebuild(link_);
		};
		
		~MarkovChain() {
			for(int site = 0; site < nSite_; ++site) {
				std::cout << site << "Up:  " << updateAcc_[2*site]/static_cast<double>(updateTot_[2*site]) << std::endl;
				std::cout << site << "Down:  " << updateAcc_[2*site + 1]/static_cast<double>(updateTot_[2*site + 1]) << std::endl;
			};
			
			std::cout << "Flip:  " << updateFlipAcc_/static_cast<double>(updateFlipTot_) << std::endl;
			
			delete bath_;	

			json jConfig;
			int key = 5345433; 
			jConfig["beta"] = beta_;
			jConfig["nSite"] = nSite_;
			jConfig["key"] = key;
			for(int site = 0; site < nSite_; ++site) {
				trace_[site]->saveConfig(jConfig,site);
				delete trace_[site];
			};

			IO::writeJsonToFile(simulation_.outputFolder() + "config_" + std::to_string(node_) + ".json",jConfig,false);		
		};
	private:
		Ut::Simulation& simulation_;
		int const node_;
		Ut::EngineType rng_;
		Ut::UniformDistribution ud_;
		Ut::UniformRngType urng_;
		
		double const beta_;
		double const probFlip_;
		int const nSite_;
		
		Link::Link link_;
		std::vector<Tr::Trace*> trace_;
		Ba::Bath* bath_;
		
		int signTrace_;
		int signBath_;
		
		double accSign_;
		std::valarray<double> pK_;
		Tr::Meas acc_;
		std::valarray<Ut::complex> accChiij_;//This valarray has a weird structure. It contains Chiij for all bosonic Matsubara frequencies. Sorted first by sites indices and then Matsubara frequencies.
											//accChiij_[(i*nSite_ + j)*acc_.chi.size() + n] give Chi(i,j,i\oemga_n)
		std::valarray<double> accChiijReal_;//As we can only save real numbers, we need a second array to store it
		
		std::vector<double> updateAcc_;
		std::vector<double> updateTot_;
		int updateFlipAcc_;
		int updateFlipTot_;
			
		
		void insert(int site, int spin) {
			Tr::Trace& trace = *trace_[site];
			Ba::Bath& bath = *bath_;
				
			double prob = .0;
			prob += trace.insert(spin, urng_);
			prob += bath.insert(site, spin, trace.op().time(), trace.op().ptr(), trace.opDagg().time(), trace.opDagg().ptr(), link_);

			if(std::log(1. - urng_()) < prob) { 
				signTrace_ *= trace_[site]->acceptInsert(spin);
			    signBath_ *= bath.acceptInsert();

				++updateAcc_[2*site + spin];
				
				return;
			}
			
			trace.rejectInsert(spin);
		};
		
		void erase(int site, int spin) {
			Tr::Trace& trace = *trace_[site];
			Ba::Bath& bath = *bath_;

			double prob = .0;
			prob += trace.erase(spin, urng_);			
			prob += bath.erase(site, spin, trace.op().time(), trace.opDagg().time());
			
			if(std::log(1. - urng_()) < prob) { 
				signTrace_ *= trace_[site]->acceptErase(spin);
				signBath_ *= bath.acceptErase();

				++updateAcc_[2*site + spin];
							
				return;
			}
			
			trace.rejectErase(spin);
		};
		
		void flipSpin(int site) {
		    trace_[site]->flip();
		    if(!tryFlip()) trace_[site]->flip();
		};
		
		void flipSite(int siteA, int siteB) {
			//Don't allow the same site (siteB is in [0,nSite_-2])
			if(siteA == siteB){
				siteB = nSite_-1;
			}
			std::swap(trace_[siteA], trace_[siteB]);
			if(!tryFlip()) std::swap(trace_[siteA], trace_[siteB]);
		};
		
		int tryFlip() {
			Ba::Bath* newBath = new Ba::Bath();
			
			for(int site = 0; site < nSite_; ++site) {
				Tr::Trace& trace = *trace_[site];
				
				for(int spin = 0; spin < 2; ++spin)
				for(Tr::Trace::Operators::const_iterator it = trace.operators(spin).begin(); it != trace.operators(spin).end(); ++it) 
					it->type() ? newBath->addDagg(site, spin, it->time(), it->ptr()) : newBath->add(site, spin, it->time(), it->ptr());
			}
		
			newBath->rebuild(link_);
		
			if(urng_() < std::abs(newBath->det()/bath_->det())) {
				signTrace_ = 1;
				for(std::vector<Tr::Trace*>::iterator it = trace_.begin(); it != trace_.end(); ++it)
					signTrace_ *= (*it)->sign();
				
				delete bath_; bath_ = newBath;
				
				signBath_ = bath_->det() > .0 ? 1 : -1;
				
				++updateFlipAcc_;
				
				return 1;
			};
			delete newBath;
			
			return 0;
		};
	};
};

#endif




















