#ifndef __NEWIO
#define __NEWIO

#include <json_spirit.h>
#include <fstream>
#include "Patrick/Utilities.h"
#include "Patrick/Hyb.h"
namespace newIO
{

	//Used as a generic type that all type are children of
	struct GenericReader{
		GenericReader(json_spirit::Value_type valueType) : valueType_(valueType){};
		virtual void write(json_spirit::mObject&,std::string const&) const = 0;
		virtual ~GenericReader() = default;
		json_spirit::Value_type getType(){
			return valueType_;
		}
		std::string getTypeName() const{
			switch(valueType_){
					case json_spirit::str_type: return "String"; break;
					case json_spirit::real_type: return "Double"; break;
					case json_spirit::int_type: return "Integer"; break;
					case json_spirit::array_type: return "ArrayDouble"; break;
					case json_spirit::obj_type: return "Green"; break;
					default : return "Unsupported type"; break;
				}
		}
		virtual std::string getString() const
		{
			throw std::runtime_error("Calling getString on a " + getTypeName() + "reader in newIO.h");
			return "";
		}
		virtual int getInt() const
		{
			throw std::runtime_error("Calling getInt on a " + getTypeName() + "reader in newIO.h");
			return 0;
		}
		virtual double getDouble(int k) const
		{
			throw std::runtime_error("Calling getDouble(k) on a " + getTypeName() + "reader in newIO.h");
			return 0.0;
		}
		virtual double getDouble() const
		{
			throw std::runtime_error("Calling getDouble() on a " + getTypeName() + "reader in newIO.h");
			return 0.0;
		}
		virtual std::size_t getSize() const {
			throw std::runtime_error("Calling getSize on a " + getTypeName() + "reader in newIO.h");
			return 0;
		}
		virtual double getFM() const {
			throw std::runtime_error("Calling getFM on a " + getTypeName() + "reader in newIO.h");
			return 0.0;
		}
		virtual double getSM() const {
			throw std::runtime_error("Calling getSM on a " + getTypeName() + "reader in newIO.h");
			return 0.0;
		}
		virtual double getBeta() const {
			throw std::runtime_error("Calling getBeta on a " + getTypeName() + "reader in newIO.h");
			return 0.0;
		}
		virtual std::complex<double> getFunction(int k) const {
			throw std::runtime_error("Calling getFunction on a " + getTypeName() + "reader in newIO.h");
			return std::complex<double>(0.0,0.0);
		}


		virtual void setDouble(double val){
			throw std::runtime_error("Calling setDouble on a " + getTypeName() + "reader in newIO.h");
		};
		virtual void setString(std::string val){
			throw std::runtime_error("Calling setString on a " + getTypeName() + "reader in newIO.h");
			
		};
		virtual void addSign(double sign){
			throw std::runtime_error("Oups, adding the sign is not implemented yet on a " + getTypeName() + " type object in newIO.h");
		}





	private:
		const json_spirit::Value_type valueType_;

	};

	struct StringReader : public GenericReader{
		StringReader(json_spirit::mValue const& jEntry,json_spirit::Value_type valueType): GenericReader(valueType), string_(jEntry.get_str()) {}; 
		std::string getString() const {
			return string_;
		};
		void setString(std::string val){
			string_ = val;
		};
		void write(json_spirit::mObject& jHyb,std::string const& entryName) const
		{
			jHyb[entryName] = string_;
		}
	private:
		std::string string_;
	};

	struct IntReader : public GenericReader{
		IntReader(json_spirit::mValue const& jEntry,json_spirit::Value_type valueType): GenericReader(valueType),value_(jEntry.get_int()) {}; 
		IntReader(int value) : GenericReader(json_spirit::int_type), value_(value) {}; 
		int getInt() const {
			return value_;
		};
		double getDouble() const {
			return static_cast<double>(value_);
		};
		void setInt(int val){
			value_ = val;
		};
		void write(json_spirit::mObject& jHyb,std::string const& entryName) const
		{
			jHyb[entryName] = value_;
		}
		void addSign(double sign){
			value_/=sign;
		}
	private:
		int value_;
	};

	struct DoubleReader : public GenericReader{
		DoubleReader(json_spirit::mValue const& jEntry,json_spirit::Value_type valueType): GenericReader(valueType),value_(jEntry.get_real()) {}; 
		DoubleReader(double value) : GenericReader(json_spirit::real_type),value_(value) {}; 
		double getDouble() const {
			return value_;
		};
		int getInt() const {
			return static_cast<int>(value_);
		};
		void setDouble(double val){
			value_ = val;
		};
		void write(json_spirit::mObject& jHyb,std::string const& entryName) const
		{
			jHyb[entryName] = value_;
		}
		void addSign(double sign){
			value_/=sign;
		}
	private:
		double value_;
	};

	struct MeanErrorReader : public GenericReader{
		MeanErrorReader(json_spirit::mValue const& jEntry,json_spirit::Value_type valueType): GenericReader(valueType),mean_(jEntry.get_obj().at("mean").get_real()), error_(jEntry.get_obj().at("error").get_real()) {}; 
		double getMean() const {
			return mean_;
		};
		double getError() const{
			return error_;
		};
		void write(json_spirit::mObject& jHyb,std::string const& entryName) const
		{
			json_spirit::mObject jObject; 
			jObject["mean"] = mean_;
			jObject["error"] = error_;
			jHyb[entryName] = jObject;
		}
	private:
		double mean_;
		double error_;
	};

	struct ArrayDoubleReader : public GenericReader{
		ArrayDoubleReader(json_spirit::mValue const& jEntry,json_spirit::Value_type valueType): GenericReader(valueType){
			json_spirit::mArray const& jArray = jEntry.get_array();
			for(std::size_t i=0;i<jArray.size();i++)
			{
				values_.push_back(newIO::DoubleReader(jArray[i].get_real()));
			}
		}; 
		double getDouble(int k) const {
			return values_[k].getDouble();
		};
		std::size_t getSize() const {
			return values_.size();
		};
		void write(json_spirit::mObject& jHyb,std::string const& entryName) const
		{
			json_spirit::mArray jArray; 
			for(std::size_t i=0;i<values_.size();i++)
			{
				jArray.push_back(getDouble(i));
			}
			jHyb[entryName] = jArray;
		}
		void addSign(double sign){
			for(std::vector<DoubleReader>::iterator it = values_.begin();it != values_.end();it++)
			{
				it->addSign(sign);
			}
		}
	private:
		std::vector<DoubleReader> values_;
	};

	struct GreenReader : public GenericReader{
		GreenReader(json_spirit::mValue const& jEntry,json_spirit::Value_type valueType): GenericReader(valueType){
			json_spirit::mObject const& jObject = jEntry.get_obj();
			FM_ = jObject.at("First Moment").get_real();
			SM_ = jObject.at("Second Moment").get_real();
			beta_ = jObject.at("beta").get_real();
			Hyb::read(jObject.at("real").get_array(),jObject.at("imag").get_array(), beta_, values_, beta_);
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
		void write(json_spirit::mObject& jHyb,std::string const& entryName) const
		{
			json_spirit::mObject jObject;
			jObject["First Moment"] = FM_;
			jObject["Second Moment"] = SM_;
			jObject["beta"] = beta_;
			json_spirit::mArray real; 
			json_spirit::mArray imag; 
			for(std::size_t i=0;i<values_.size();i++)
			{
				std::complex<double> function = getFunction(i);
				real.push_back(function.real());
				imag.push_back(function.imag());
			}
			jObject["real"] = real;
			jObject["imag"] = imag;
			jHyb[entryName] = jObject;
		}
		void addSign(double sign){
			for(std::vector<std::complex<double> >::iterator it = values_.begin();it != values_.end();it++)
			{
				*it /=sign;
			}
		}
	private:
		double FM_;
		double SM_;
		double beta_;
		std::vector<std::complex<double> > values_;
	};

	//We copy IO.h's function to be able to differenciate between types and get every type of read right
	struct GenericReadFunc {
		GenericReadFunc(std::string fileName,std::string nodeName){
			json_spirit::mObject jHybComplete;
			std::ifstream file(fileName.c_str()); 
			if(file) {
				json_spirit::mValue temp;
				json_spirit::read(file, temp); 
				jHybComplete = temp.get_obj();
			}else{
				throw std::runtime_error(fileName + " not found.");
			}
			json_spirit::mObject jHyb;
			if(!nodeName.empty())
			{
				jHyb = jHybComplete.at(nodeName).get_obj();
			}else
			{
				jHyb = jHybComplete;
			}
			int index = 0;
			std::cout << "Reading from " << fileName << " at " << nodeName << std::endl;
			for(json_spirit::mObject::const_iterator it = jHyb.begin(); it != jHyb.end(); ++it, ++index) {
				index_[it->first] = index;
				//Now we have to create a custom object for each type of array
				json_spirit::mValue const& jEntry = it->second;
				json_spirit::Value_type const& entryType = jEntry.type();
				GenericReader* container = 0;
				switch(entryType){
					case json_spirit::str_type: container = new StringReader(jEntry,entryType); break;
					case json_spirit::real_type: container = new DoubleReader(jEntry,entryType); break;
					case json_spirit::int_type: container = new IntReader(jEntry,entryType); break;
					case json_spirit::array_type: 
						if(jEntry.get_array().size() == 1)
						{
							container = new DoubleReader(jEntry.get_array()[0],entryType);
						}else
						{
							container = new ArrayDoubleReader(jEntry,entryType);
						}
						break;
					case json_spirit::obj_type: container = new GreenReader(jEntry,entryType); break;
					default :break;
				}
				std::cout << it->first.c_str() << " ";
				hyb_.push_back(container);
			}
			std::cout << std::endl;
		};
		struct GenericReader* operator()(std::string str){ 
			if(index_.find(str) == index_.end())
				throw std::runtime_error("Entry " + str + " not found.");
			return hyb_[index_.at(str)];
		};
		struct GenericReader* operator()(std::string str,bool& found){ 
			if(index_.find(str) == index_.end())
			{
				return NULL;
			}
			return hyb_[index_.at(str)];
		};
		void write(std::string fileName)
		{
			json_spirit::mObject jHyb;
			for(std::map<std::string, std::size_t>::iterator it = index_.begin(); it != index_.end(); ++it) {
				hyb_[it->second]->write(jHyb,it->first);
			};

			std::ofstream paramsNext(fileName); 
			json_spirit::write(jHyb, paramsNext, json_spirit::pretty_print | json_spirit::single_line_arrays | json_spirit::remove_trailing_zeros);
			paramsNext.close();
		}





		/*
		* void addSign(double sign)
		*
		* Arguments : 	sign : Computed sign of the monte Carlo Solver, to get the right data out of the program
		*
		* Description : 
		*		This function compute the real data from the signed data that come out of the impurity solver program.
		*		It is mandatory to execute this function on the measures coming out of the Impurity Solver (Alps did it on its own and now it's not done anymore in the without alps program)
		*/
		void addSign(double sign)
		{
			for(std::map<std::string, std::size_t>::iterator it = index_.begin(); it != index_.end(); ++it) {
				hyb_[it->second]->addSign(sign);
			};
		}
		~GenericReadFunc() {
			for(std::size_t i = 0; i < index_.size(); ++i) delete hyb_[i];
		};
	private:
		std::map<std::string, std::size_t> index_;
		std::vector<GenericReader*> hyb_;
	};

}
#endif