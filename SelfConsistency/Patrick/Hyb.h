#ifndef __HYB
#define __HYB

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include "../nlohmann_json.hpp"
#include "Utilities.h"
using json=nlohmann::json;

namespace Hyb {
	/* When changing beta, the matsubara frequencies change, and the indices don't correspond to the same frequencies */
	/* This function adapts the Matsubara data contained in (jComponentIn["real"] + i*jComponentIn["imag"]) */
	/* 	from jComponentIn["beta"] to betaOut and returns the result */
	/* The result conserves the other attributes of jComponentIn */
	json read(json const& jComponentIn, double betaOut) {
		json jComponentOut = jComponentIn;
		double const betaIn = jComponentIn["beta"];
		if(jComponentIn["real"].size() != jComponentIn["imag"].size()) 
			throw std::runtime_error("Missmatch of real and imaginary array lengths.");
		
		jComponentOut["real"] = json::array();
		jComponentOut["imag"] = json::array();
		int n = 0;
		for(int m = 1; m < static_cast<int>(jComponentIn["real"].size()); ++m) { 
			double z0 = M_PI*(2*m - 1)/betaIn;
			double z1 = M_PI*(2*m + 1)/betaIn;
			
			Ut::complex h0(jComponentIn["real"][m - 1], jComponentIn["imag"][m - 1]);
			Ut::complex h1(jComponentIn["real"][m], jComponentIn["imag"][m]);
			
			for(; (2*n + 1)/betaOut <= (2*m + 1)/betaIn; ++n) {
				double w = M_PI*(2*n + 1)/betaOut;
				std::complex<double> value = h0 + (h1 - h0)*(w - z0)/(z1 - z0);
				jComponentOut["real"].push_back(value.real());
				jComponentOut["imag"].push_back(value.imag());
			}
		}
		jComponentOut["beta"] = betaOut;
		return jComponentOut;
	};
	/* Constructs the Hybridization Function in the imaginary time space */
	struct Function {
		Function(std::string name, json const& jNumericalParams, json const& jEntry) : 
		beta_(jNumericalParams["beta"]), 
		nItH_(jNumericalParams["EHyb"].get<double>()*jNumericalParams["EHyb"].get<double>()*beta_ + 1), 
		hyb_(new double[nItH_ + 1]) {
			std::cout << "Initialising data entry " << name << " ... ";
			
			double const FM = jEntry["First Moment"]; 
			double const SM = jEntry["Second Moment"];
			
			for(int i = 0; i < nItH_ + 1; ++i) {
				double time = beta_*i/static_cast<double>(nItH_);
				hyb_[i] = -FM/2. + SM*(time - beta_/2.)/2.;
			}
			   
			/* Adapts the Hyb entry to the current beta_ */
			json jHybWithRightBeta = read(jEntry, beta_);
			json jReal = jHybWithRightBeta["real"];
			json jImag = jHybWithRightBeta["imag"];
			/* Creates the imaginary time data */
			for(unsigned int m = 0; m < jReal.size(); ++m) {
				Ut::complex z(.0, M_PI*(2*m + 1)/beta_);
				Ut::complex value = std::complex<double>(jReal[m],jImag[m]) - FM/z - SM/(z*z);
				for(int i = 0; i < nItH_ + 1; ++i) {
					double arg = M_PI*(2*m + 1)*i/static_cast<double>(nItH_); 
					hyb_[i] += 2./beta_*(value.real()*std::cos(arg) + value.imag()*std::sin(arg));
				}
			}
			
			std::cout << "Ok" << std::endl;
		};
		
		double get(double time) const { 
			double it = time/beta_*nItH_; int i0 = static_cast<int>(it);
			return (1. - (it - i0))*hyb_[i0] + (it - i0)*hyb_[i0 + 1];
		};
		
		~Function() { delete[] hyb_;};
	private:
		double const beta_;
		int const nItH_;
		
		double* const hyb_;
	};
	
	//------------------------------------------------------------------------------------------------------------------------------------------
	/**************************************************************************************************/
	/*              This function adapts the Matsubara data file completely to betaOut                */
	/*                          The jObject is just changed in place                                  */
	/**************************************************************************************************/
	void accountForDifferentBeta(json& jObject,double const betaOut){
		json const& jTemp = jObject;
		/* It changes the beta for all components in jObject */
		for (auto& el : jTemp.items()){
			std::cout << "Initialising data entry " << el.key() << " ... ";
			jObject[el.key()] = read(jTemp[el.key()],betaOut);
			std::cout << "Ok" << std::endl;
		} 
	}
};

#endif
