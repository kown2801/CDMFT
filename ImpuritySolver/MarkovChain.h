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
#include <boost/random.hpp>
#include <json_spirit.h>
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
		
		MarkovChain(json_spirit::mObject const& jNumericalParams, json_spirit::mObject const& jHyb, json_spirit::mArray const& jLinkN, json_spirit::mArray const& jLinkA, Ut::Simulation& simulation) :
		simulation_(simulation),
		node_(mpi::rank()),
		rng_(jNumericalParams.at("SEED").get_int()),
		urng_(rng_, Ut::UniformDistribution(.0, 1.)),
		beta_(jNumericalParams.at("beta").get_real()),
		probFlip_(jNumericalParams.at("PROBFLIP").get_real()),
		nSite_(jLinkN.size()),
		link_(jNumericalParams, jHyb, jLinkN, jLinkA, simulation.meas()),
		trace_(nSite_, static_cast<Tr::Trace*>(0)),
		bath_(new Ba::Bath()),
		signTrace_(1),
		signBath_(1),
		accSign_(.0),
		accChiij_(.0,nSite_*nSite_),
		pK_(.0, jNumericalParams.at("EOrder").get_int()),
		acc_(jNumericalParams),
		updateAcc_(2*nSite_, .0),
		updateTot_(2*nSite_, .0),
		updateFlipAcc_(0),
		updateFlipTot_(0) {
			if(jLinkN.size() != jLinkA.size())
				throw std::runtime_error("MarkovChain: missmatch in size of LinkN and LinkA.");
			Ut::Measurements& measurements = simulation.meas();
			std::cout << jNumericalParams.at("SEED").get_int() << std::endl;
            
			json_spirit::mValue jPreviousConfigValue;


			try{
				mpi::read_json_all_processors(simulation_.outputFolder() + "config_" + boost::lexical_cast<std::string>(node_) + ".json", jPreviousConfigValue);
				json_spirit::mObject const& jPreviousConfig = jPreviousConfigValue.get_obj();
				if(jPreviousConfig.size()) {
					double beta = jPreviousConfig.at("beta").get_real();
					if(beta != beta_) throw std::runtime_error("MarkovChain: missmatch in beta");
					
					int nSite = jPreviousConfig.at("nSite").get_int();
					if(nSite != nSite_) throw std::runtime_error("MarkovChain: missmatch in site number");
					
					for(int site = 0; site < nSite_; ++site) {
						trace_[site] = new Tr::Trace(jNumericalParams, site, measurements, jPreviousConfig);
						
						signTrace_ *= trace_[site]->sign();
						
						for(int spin = 0; spin < 2; ++spin)
						    for(Tr::Trace::Operators::const_iterator it = trace_[site]->operators(spin).begin(); it != trace_[site]->operators(spin).end(); ++it) 
							    it->type() ? bath_->addDagg(site, spin, it->time(), it->ptr()) : bath_->add(site, spin, it->time(), it->ptr());
					};
					
					signBath_ *= bath_->rebuild(link_);
					
					int key = jPreviousConfig.at("key").get_int();
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
				json_spirit::mObject jPreviousConfig = json_spirit::mObject();
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
			std::vector<double> Chiij_local(4,0.);
			for(int site = 0; site < nSite_; ++site) {
				trace_[site]->measure(sign,Chiij_local[site]);
				k += (trace_[site]->operators(0).size() + trace_[site]->operators(1).size())/2;
			};
			for(int i = 0;i<nSite_;i++){
				for(int j = 0;j<nSite_;j++){
					accChiij_[i*nSite_ + j] += sign*Chiij_local[i]*Chiij_local[j]/beta_;
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

			accChiij_ /= NAlpsMeas;
			measurements["Chiij"] << accChiij_;
			accChiij_ = .0;
			
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

			json_spirit::mObject jConfig;
			int key = 5345433; 
			jConfig["beta"] = beta_;
			jConfig["nSite"] = nSite_;
			jConfig["key"] = key;
			for(int site = 0; site < nSite_; ++site) {
				trace_[site]->saveConfig(jConfig,site);
				delete trace_[site];
			};

			std::ofstream file((simulation_.outputFolder() + "config_" + boost::lexical_cast<std::string>(node_) + ".json").c_str());			
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
		Ba::Bath* bath_;
		
		int signTrace_;
		int signBath_;
		
		double accSign_;
		std::valarray<double> accChiij_;
		std::valarray<double> pK_;
		Tr::Meas acc_;
		
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




















