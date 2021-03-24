#ifndef __LINK
#define __LINK

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "nlohmann_json.hpp"
#include "Utilities.h"
#include "Green.h"
#include "Hyb.h"

namespace Link {
	struct Link {
	/** 
		* 
		* Link(json const& jNumericalParams, json const& jHyb, json const& jLink, Ut::Measurements& measurements)
		* 
		* Parameters :	jNumericalParams : storage of all the numerical parameters of the simulation
		*				jHyb : storage of the hybridization informations (read more in README.MD)
		*				jLink : storage of the link informations (read more in README.MD)
		*				simulation : Object used to store the parameters and the measurements throughout the simulation
		* 
		* Prints :	The status of the initialization and whether it was able to read in different files
		*
		* Description: 
		*   Constructs the hyb function from the Link and Hyb files (making it possible to read the hybridation function between each sites)
		*	Reserves some space for the green function to be saved in
		* 
		*/
		Link(json const& jNumericalParams, json const& jHyb, json const& jLink, Ut::Measurements& measurements) : 
		beta_(jNumericalParams["beta"]), 
		//The Link object includes spins up and down so the number of site is half the Link array size
		nSite_(jLink.size()/2),
		greenEntry_(new Entry[4*nSite_*nSite_]), 
		hybEntry_(new Entry[4*nSite_*nSite_]), 
		multiplicity_(jHyb.size(), 0) {	
			std::map<std::string, int> entryIndex; int index = 0;
			for (auto& el : jHyb.items()){
			    entryIndex[el.key()] = index++; 
			}
				
			
			std::cout << "Reading in link file ... " << std::endl;			
			
			for(unsigned int i = 0; i < jLink.size(); ++i) {
				if(jLink[i].size() != 2*nSite_) 
					throw(std::runtime_error("Wrong row size."));
				
				for(unsigned int j = 0; j < jLink[i].size(); ++j) {
					std::string entry = jLink[i][j]; std::cout << entry << " ";						
					Entry temp;
					
					if(entry.compare("empty")) {
						if(entryIndex.find(entry) == entryIndex.end())
							throw(std::runtime_error("Entry in hybridisation not found."));
					
						
						multiplicity_[entryIndex.at(entry)] += 1;  
						temp.index = entryIndex.at(entry);
					}	
					
					temp.fact = 1.;
					temp.arg = 1.;
					//We need to take the opposite of the conjugate for the down spins as per the Nambu representation
					if (i >= nSite_ && j>= nSite_){
						temp.fact *= -1.;
						temp.arg *= -1.;
					}
					
					hybEntry_[i  + 2*nSite_*j] = greenEntry_[j  + 2*nSite_*i] = temp;
				}
			}
			
			std::cout << std::endl << "... Ok" << std::endl << std::endl;
			
			hyb_ = allocHyb_.allocate(multiplicity_.size());
			green_ = allocGreen_.allocate(multiplicity_.size());
			
			for(std::map<std::string, int>::const_iterator it = entryIndex.begin(); it != entryIndex.end(); ++it) {
				new(hyb_ + it->second) Hyb::Function(it->first, jNumericalParams, jHyb[it->first]);
				new(green_ + it->second) Green::Meas(it->first, jNumericalParams, measurements);
			}
		};
		/** 
		* 
		* double operator()(Op const& opL, Op const& opR) const
		* 
		* Parameters :	opL,opR : operators onwhich we eant the Hybridiation function
		* 
		* Return Value : The value of the hybridation function between the two operators at time time
		* 
		*/
		template<class Op> double operator()(Op const& opL, Op const& opR) const { 
			Entry entry = hybEntry_[(opL.spin()*nSite_ + opL.site()) + 2*nSite_*(opR.spin()*nSite_ + opR.site())];
			
			if(entry.index != -1) {
				double time = entry.arg*(opL.time() - opR.time());
				
				if(time < .0) {
					time += beta_;
					entry.fact *= -1;
				}
				
				return entry.fact*hyb_[entry.index].get(time);
			}
			
			return .0;
		};
				/** 
		* 
		* template<class GreenIterator>
		* void measure(int sign, GreenIterator begin, GreenIterator end)
		* 
		* Parameters :	sign : current sign of the bath
		*				begin : start of the green matrix
		*				end : end of the green matrix
		* 
		* Description :
		*	Adds the current green function value to the green function measurements
		*	For each element in the Green matrix (between begin and end), adds its contribution to the imaginary time green's function.
		* 
		*/	
		template<class GreenIterator>
		void measure(int sign, GreenIterator begin, GreenIterator end) {
			for(GreenIterator it = begin; it != end; ++it) {
				Entry entry = greenEntry_[(it.opR().spin()*nSite_ + it.opR().site()) + 2*nSite_*(it.opL().spin()*nSite_ + it.opL().site())];
				
				if(entry.index != -1) {
					double time = entry.arg*(it.opR().time() - it.opL().time()); 
					double value = sign*entry.fact*it.value();
					
					if(time < .0) {
						value *= -1.;
						time += beta_;
					}
					
					green_[entry.index].add(time, value);  //Spin !!!!
				}
			}
		};
		/** 
		* 
		* void store(Ut::Measurements& measurements, int measurementsFromLastStore)
		* 
		* Parameters :	measurements : variables used to store the measurements, used for output
		*				measurementsFromLastStore : Number of measurements done simce the last time we stored some measurements
		* 
		* Description :
		*	Stores the measurements for every green's component
		* 
		*/
		void store(Ut::Measurements& measurements, int measurementsFromLastStore) {
			for(unsigned int i = 0; i < multiplicity_.size(); ++i) 
				green_[i].measure(measurements, multiplicity_[i]*measurementsFromLastStore);
		};
		
		~Link() { 
			for(unsigned int i = 0; i < multiplicity_.size(); ++i) {
				green_[i].~Meas();
				hyb_[i].~Function();				
			}
			allocGreen_.deallocate(green_, multiplicity_.size());
		    allocHyb_.deallocate(hyb_, multiplicity_.size());	
			
			delete[] hybEntry_;
			delete[] greenEntry_; 
		};
	private:
		struct Entry { 
			Entry() : index(-1) {};
			double fact; double arg; int index;
		};
		
		std::allocator<Hyb::Function> allocHyb_;
		std::allocator<Green::Meas> allocGreen_;
		
		double const beta_;
		std::size_t const nSite_;
		Entry* const greenEntry_;
		Entry* const hybEntry_;
		
		std::vector<int> multiplicity_;
		Hyb::Function* hyb_;
		Green::Meas* green_;
	};
};

#endif