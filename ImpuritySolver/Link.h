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
		nSite_(jLink.size()),
		greenEntry_(new int[nSite_*nSite_]), 
		hybEntry_(new int[nSite_*nSite_]), 
		multiplicity_(jHyb.size(), 0) {	
			std::map<std::string, int> entryIndex; int index = 0;
			for (auto& el : jHyb.items()){
			    entryIndex[el.key()] = index++; 
			}
			
			std::cout << "Reading in link file ... " << std::endl;			
			
			for(unsigned int i = 0; i < jLink.size(); ++i) {
				if(jLink[i].size() != nSite_) 
					throw(std::runtime_error("Wrong row size."));
				
				for(unsigned int j = 0; j < jLink[i].size(); ++j) {
					std::string entry = jLink[i][j]; std::cout << entry << " ";						
					
					if(entryIndex.find(entry) == entryIndex.end())
						throw(std::runtime_error("Entry in hybridisation not found."));
					
					multiplicity_[entryIndex.at(entry)] += 2;    //Spin !!!!!
					
					hybEntry_[i  + nSite_*j] = greenEntry_[j  + nSite_*i] = entryIndex.at(entry);						
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
		* double operator()(int siteDagg, int site, double time) const
		* 
		* Parameters :	siteDagg : index of the site where a particle is created
		*				site : index of the site where a particle is annihilated
		*				time : time at which the hybridizatio function should be returnes
		* 
		* Return Value : The value of the hybridation function between the two sites at time time
		* 
		*/
		double operator()(int siteDagg, int site, double time) const { 
			int sign = 1;
			
			if(time < .0) {
				time += beta_;
				sign = -1;
			}
			
			return sign*hyb_[hybEntry_[siteDagg + nSite_*site]].get(time);
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
				double time = it.op().time() - it.opDagg().time(); 
				double value = sign*it.value();
				
				if(time < .0) {
					value *= -1.;
					time += beta_;
				}
				
				green_[greenEntry_[it.op().site() + nSite_*it.opDagg().site()]].add(time, value);  //Spin !!!!
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
				green_[i].store(measurements, multiplicity_[i]*measurementsFromLastStore);
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
		std::allocator<Hyb::Function> allocHyb_;
		std::allocator<Green::Meas> allocGreen_;
		
		double const beta_;
		std::size_t const nSite_;
		int* const greenEntry_;
		int* const hybEntry_;
		
		std::vector<int> multiplicity_;
		Hyb::Function* hyb_;
		Green::Meas* green_;
	};
};

#endif