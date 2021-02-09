#ifndef __IO
#define __IO

#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include "Utilities.h"
#include "Hyb.h"
#include "../newIO.h"

namespace IO {
		
	template<class F>
	struct GenericWriteFunc {
		GenericWriteFunc() {};
		F& operator()(std::string str) { return hyb_[str];};
		std::map<std::string, F>& getMap(){return hyb_;};
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
	
	typedef GenericReadFunc<Hyb::Read> ReadFunc;
	typedef GenericWriteFunc<Hyb::Write> WriteFunc;
	typedef GenericWriteFunc<newIO::DatWriter> WriteDat;
};

#endif