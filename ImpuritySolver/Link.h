#ifndef __LINK
#define __LINK

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <json_spirit.h>
#include "Utilities.h"
#include "Green.h"
#include "Hyb.h"

namespace Link {
	struct Link {
		Link(json_spirit::mObject const& jNumericalParams, json_spirit::mObject const& jHyb, json_spirit::mArray const& jLinkN, json_spirit::mArray const& jLinkA, Ut::Measurements& measurements) : 
		beta_(jNumericalParams.at("beta").get_real()), 
		nSite_(jLinkN.size()),
		greenEntry_(new Entry[4*nSite_*nSite_]), 
		hybEntry_(new Entry[4*nSite_*nSite_]), 
		multiplicity_(jHyb.size(), 0) {	
			std::map<std::string, int> entryIndex; int index = 0;
			for(json_spirit::mObject::const_iterator it = jHyb.begin(); it != jHyb.end(); ++it) 
				entryIndex[it->first] = index++; 
			
			std::cout << "Reading in link file ... " << std::endl;			
			
			for(unsigned int i = 0; i < jLinkN.size(); ++i) {
				json_spirit::mArray jRow = jLinkN[i].get_array();
				if(jRow.size() != nSite_) 
					throw(std::runtime_error("Wrong row size."));
				
				for(unsigned int j = 0; j < jRow.size(); ++j) {
					std::string entry = jRow[j].get_str(); std::cout << entry << " ";						
					
					if(entryIndex.find(entry) == entryIndex.end())
						throw(std::runtime_error("Entry in hybridisation not found."));
					
					Entry temp;
					
					multiplicity_[entryIndex.at(entry)] += 2;    //Spin !!!!!
					temp.index = entryIndex.at(entry);
					
					temp.fact = 1.;
					temp.arg = 1.;
					
					hybEntry_[i  + 2*nSite_*j] = greenEntry_[j  + 2*nSite_*i] = temp;
					
					temp.fact *= -1.;
					temp.arg *= -1.;
					
					hybEntry_[(i + nSite_)  + 2*nSite_*(j + nSite_)] = greenEntry_[(j + nSite_)  + 2*nSite_*(i + nSite_)] = temp;
				}
			}
			
			for(unsigned int i = 0; i < jLinkA.size(); ++i) {
				json_spirit::mArray jRow = jLinkA[i].get_array();
				if(jRow.size() != nSite_) 
					throw(std::runtime_error("Wrong row size."));
				
				for(unsigned int j = 0; j < jRow.size(); ++j) {
					std::string entry = jRow[j].get_str(); std::cout << entry << " ";						
					
					Entry temp;
					temp.fact = 1.;
					temp.arg = 1.;
					
					if(entry.compare("empty")) {
						if(entryIndex.find(entry) == entryIndex.end())
							throw(std::runtime_error("Entry in hybridisation not found."));
						
						multiplicity_[entryIndex.at(entry)] += 2;    //Spin !!!!!
						temp.index = entryIndex.at(entry);
					} 
					
					hybEntry_[(i + nSite_)  + 2*nSite_*j] = greenEntry_[(j + nSite_)  + 2*nSite_*i] = temp;
					hybEntry_[i + 2*nSite_*(j + nSite_)] = greenEntry_[j  + 2*nSite_*(i + nSite_)] = temp;
				}
			}
			
			std::cout << std::endl << "... Ok" << std::endl << std::endl;
			
			hyb_ = allocHyb_.allocate(multiplicity_.size());
			green_ = allocGreen_.allocate(multiplicity_.size());
			
			for(std::map<std::string, int>::const_iterator it = entryIndex.begin(); it != entryIndex.end(); ++it) {
				new(hyb_ + it->second) Hyb::Function(it->first, jNumericalParams, jHyb.at(it->first).get_obj());
				new(green_ + it->second) Green::Meas(it->first, jNumericalParams, measurements);
			}
		};
		
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
		
		void measure(Ut::Measurements& measurements, int NAlpsMeas) {
			for(unsigned int i = 0; i < multiplicity_.size(); ++i) 
				green_[i].measure(measurements, multiplicity_[i]*NAlpsMeas);
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