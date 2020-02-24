#ifndef __NEWIO
#define __NEWIO

#include <json_spirit.h>
#include <fstream>
#include "Patrick/Utilities.h"
#include "Patrick/Hyb.h"
namespace newIO
{

	//Used as a generic type that all type are children of
	struct GenericReader{};

	struct StringReader : public GenericReader{
		StringReader(json_spirit::mObject const& jEntry): string_(jEntry.begin()->second.get_str()) {}; 
		const std::string& getString() const {
			return string_;
		};
		void setString(std::string val){
			string_ = val;
		};
	private:
		std::string string_;
	};

	struct DoubleReader : public GenericReader{
		DoubleReader(json_spirit::mObject const& jEntry) : value_(jEntry.begin()->second.get_real()) {}; 
		double getDouble() const {
			return value_;
		};
		void setDouble(double val){
			value_ = val;
		};
	private:
		double value_;
	};

	struct MeanErrorReader : public GenericReader{
		MeanErrorReader(json_spirit::mObject const& jEntry) : mean_(jEntry.at("mean").get_real()), error_(jEntry.at("error").get_real()) {}; 
		double getMean() const {
			return mean_;
		};
		double getError() const{
			return error_;
		};
	private:
		double mean_;
		double error_;
	};

	struct ArrayDoubleReader : public GenericReader{
		ArrayDoubleReader(json_spirit::mObject const& jEntry){
			json_spirit::mArray const& jArray = jEntry.begin()->second.get_array();
			for(std::size_t i=0;i<jEntry.size();i++)
			{
				values_.push_back(newIO::MeanErrorReader(jArray[i].get_obj()));
			}
		}; 
		double getMean(int k) const {
			return values_[k].getMean();
		};
		double getError(int k) const  {
			return values_[k].getError();
		};
		std::size_t getSize() const {
			return values_.size();
		};
	private:
		std::vector<MeanErrorReader> values_;
	};

	struct GreenReader : public GenericReader{
		GreenReader(json_spirit::mObject const& jEntry): FM_(jEntry.at("First Moment").get_real()), SM_(jEntry.at("Second Moment").get_real()),beta_(jEntry.at("beta").get_real()){
				Hyb::read(jEntry.at("real").get_array(),jEntry.at("imag").get_array(), beta_, values_, beta_);
		};
		double getFM() const {
			return FM_;
		};
		double getSM() const {
			return SM_;
		};
		double getBeta() const {
			return beta_;
		};
		std::complex<double> getFunction(int k) const {
			return values_[k];
		};
		std::size_t getSize() const {
			return values_.size();
		};
	private:
		double FM_;
		double SM_;
		double beta_;
		std::vector<std::complex<double> > values_;
	};

	//We copy IO.h's function to be able to differenciate between types and get every type of read right
	struct GenericReadFunc {
		GenericReadFunc(std::string fileName){
			json_spirit::mObject jHyb;
			std::ifstream file(fileName.c_str()); 
			if(file) {
				json_spirit::mValue temp;
				json_spirit::read(file, temp); 
				jHyb = temp.get_obj();
			} else 
				throw std::runtime_error(fileName + " not found.");
			
			int index = 0;
			for(json_spirit::mObject::const_iterator it = jHyb.begin(); it != jHyb.end(); ++it, ++index) {
				index_[it->first] = index;
				//Now we have to create a custom object for each type of array
				json_spirit::mObject const& jEntry = it->second.get_obj();
				json_spirit::Value_type const& firstEntryType = jEntry.begin()->second.type();
				GenericReader* container;
				switch(firstEntryType){
					case json_spirit::str_type: container = new StringReader(jEntry); break;
					case json_spirit::real_type: container = new DoubleReader(jEntry); break;
					case json_spirit::array_type: container = new ArrayDoubleReader(jEntry); break;
					case json_spirit::obj_type:
						if(!(jEntry.at("mean").is_null()))
						{
							//The mean error case 
							container = new MeanErrorReader(jEntry);
						}else
						{
							//The Green's fucntion case
							container = new GreenReader(jEntry);

						}
					break;
					default :return;
				}
				std::cout << it->first << std::endl;
				hyb_.push_back(container);
			}
		};
		struct GenericReader* operator()(std::string str) { 
			if(index_.find(str) == index_.end())
				throw std::runtime_error("Entry " + str + " not found.");
			return hyb_[index_.at(str)];
		};
		struct GenericReader* operator()(std::string str,bool& found){ 
			if(index_.find(str) == index_.end())
			{
				found = false;
				return hyb_[index_.begin()->second];
			}
			found = true;
			return hyb_[index_.at(str)];
		};
		~GenericReadFunc() {
			for(std::size_t i = 0; i < index_.size(); ++i) delete hyb_[i];
		};
	private:
		std::map<std::string, std::size_t> index_;
		std::vector<GenericReader*> hyb_;
	};
}
#endif