#ifndef __IO
#define __IO

#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <alps/alea.h>
#include <alps/scheduler/montecarlo.h>
#include <json_spirit.h>
#include "Utilities.h"
#include "Hyb.h"
#include "Green.h"

namespace IO {
	struct Scal {
		typedef Ut::MeasEntry<double> RET;
		
		Scal() : isEntry_(false) {};
		RET const& operator()(std::string name, alps::scheduler::MCSimulation& sim) {
			if(!isEntry_) {
			    alps::RealObsevaluator temp = sim.get_measurements()[name];
				entry_ = Ut::MeasEntry<double>(temp.mean(), temp.error());
				isEntry_ = true;
			}
			
			return entry_;
		}
	private:
		bool isEntry_; 
		RET entry_; 
	};
	
	struct Vec {
		typedef std::vector<Ut::MeasEntry<double> > RET;
		
		Vec() {};
		RET const& operator()(std::string name, alps::scheduler::MCSimulation& sim) {	
			if(!vec_.size()) {
				alps::RealVectorObsevaluator temp = sim.get_measurements()[name];
				
				std::valarray<double> mean = temp.mean(); 
				std::valarray<double> error = temp.error();
				
				for(std::size_t i = 0; i < mean.size(); ++i)  
					vec_.push_back(Ut::MeasEntry<double>(mean[i], error[i]));
			}
			
			return vec_;
        }
	private:
		RET vec_;
	};
	
	template<class O>
	struct GenericReadObs {
		GenericReadObs(alps::scheduler::MCSimulation& sim) : sim_(sim) {};
		typename O::RET const& operator()(std::string str) {
			return o_[str](str, sim_);
		};
	private:
		alps::scheduler::MCSimulation& sim_;
		std::map<std::string, O> o_;
	};
	
	template<class F>
	struct GenericWriteFunc {
		GenericWriteFunc() {};
		F& operator()(std::string str) { return hyb_[str];};
		void write(double beta, std::string fileName) {
			json_spirit::mObject jHyb;
			for(typename std::map<std::string, F>::iterator it = hyb_.begin(); it != hyb_.end(); ++it) {
				json_spirit::mObject temp; 
				it->second.write(beta, temp);
				jHyb[it->first] = temp;
			};
			std::ofstream file(fileName.c_str());			
			json_spirit::write(jHyb, file, json_spirit::pretty_print | json_spirit::single_line_arrays);			
			file.close();
		};
	private:
		std::map<std::string, F> hyb_;
	};
	
	template<class F>
	struct GenericReadFunc {
		GenericReadFunc(std::string fileName) {
			json_spirit::mObject jHyb;
			std::ifstream file(fileName.c_str()); 
			if(file) {
				json_spirit::mValue temp;
				json_spirit::read(file, temp); 
				jHyb = temp.get_obj();
			} else 
				throw std::runtime_error(fileName + " not found.");
			
			int index = 0;
			hyb_ = alloc_.allocate(jHyb.size());
			for(json_spirit::mObject::const_iterator it = jHyb.begin(); it != jHyb.end(); ++it, ++index) {
				index_[it->first] = index;
				new(hyb_ + index) F(it->second.get_obj());
			};
		};
		F const& operator()(std::string str) { 
			if(index_.find(str) == index_.end())
				throw std::runtime_error("Entry " + str + " not found.");
			return hyb_[index_.at(str)];
		};
		~GenericReadFunc() {
			for(std::size_t i = 0; i < index_.size(); ++i) hyb_[i].~F();
			alloc_.deallocate(hyb_, index_.size());
		};
	private:
		std::allocator<F> alloc_;
		std::map<std::string, std::size_t> index_;
		F* hyb_;
	};
	
	//----------------------------------------------------------------------------------------------------------------------------------------------------------
	
	typedef GenericReadObs<Scal> ReadScal;
	typedef GenericReadObs<Vec> ReadVec;
	typedef GenericReadObs<Green::Entry> ReadGreen;
	
	typedef GenericReadFunc<Hyb::Read> ReadFunc;
	typedef GenericWriteFunc<Hyb::Write> WriteFunc;
};

#endif