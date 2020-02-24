#ifndef __HYB
#define __HYB

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <json_spirit.h>
#include "Utilities.h"

namespace Hyb {
	void read(json_spirit::mArray const& jReal, json_spirit::mArray const& jImag, double betaIn, std::vector<Ut::complex>& Out, double betaOut) {
		if(jReal.size() != jImag.size()) 
			throw std::runtime_error("Missmatch of real and imaginary array lengths.");
		
		int n = 0;
		for(int m = 1; m < static_cast<int>(jReal.size()); ++m) { 
			double z0 = M_PI*(2*m - 1)/betaIn;
			double z1 = M_PI*(2*m + 1)/betaIn;
			
			Ut::complex h0(jReal[m - 1].get_real(), jImag[m - 1].get_real());
			Ut::complex h1(jReal[m].get_real(), jImag[m].get_real());
			
			for(; (2*n + 1)/betaOut <= (2*m + 1)/betaIn; ++n) {
				double w = M_PI*(2*n + 1)/betaOut;
				Out.push_back(h0 + (h1 - h0)*(w - z0)/(z1 - z0));
			}
		}
	};
	
	
	//------------------------------------------------------------------------------------------------------------------------------------------
	
	struct Read : public std::vector<Ut::complex> {
		Read(json_spirit::mObject const& jEntry) {
			std::cout << "Reading in data file ... ";

			beta_ = jEntry.at("beta").get_real(); 
			FM_ = jEntry.at("First Moment").get_real();
			SM_ = jEntry.at("Second Moment").get_real();
			
			read(jEntry.at("real").get_array(), jEntry.at("imag").get_array(), jEntry.at("beta").get_real(), *this, beta_);
			
			std::cout << "Ok" << std::endl;
		};
			
		double beta() const { return beta_;};
		double FM() const { return FM_;};
		double SM() const { return SM_;};				
	private:
		double beta_;
		double FM_;
		double SM_;
	};
	
	//---------------------------------------------------------------------------------------------------------------------------------
	
	struct Write : public std::vector<Ut::complex> {
		Write() : FM_(.0), SM_(.0) {};
		double& FM() { return FM_;};
		double& SM() { return SM_;};
		
		void write(double beta, json_spirit::mObject& jEntry) {	
			std::cout << "Writing function ... ";
			
			jEntry["beta"] = beta; 
			jEntry["First Moment"] = FM_;  
			jEntry["Second Moment"] = SM_;
            			
			std::vector<double> real(size());
			std::vector<double> imag(size());			
			for(std::size_t n = 0; n < size(); ++n) {
				real[n] = at(n).real();
				imag[n] = at(n).imag();
 			}
			
			jEntry["real"] = json_spirit::mValue(real.begin(), real.end()); 
			jEntry["imag"] = json_spirit::mValue(imag.begin(), imag.end()); 
	
			std::cout << "Ok" << std::endl;
		};
	private:
		double FM_;
		double SM_;
	};
};

#endif