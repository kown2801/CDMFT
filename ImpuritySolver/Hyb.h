#ifndef __HYB
#define __HYB

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include "nlohmann_json.hpp"
#include "Utilities.h"
using json=nlohmann::json;

namespace Hyb {

	/** 
	* 
	* json read(json const& jComponentIn, double betaOut)
	* 
	* Parameters :	jComponentIn : component of the Matsubara function that is being read
	*				betaOut : beta of the program in which the function needs to be transformed
	*
	* Returns:
	*	A json object containing the real part ("real" attribute) and the imaginary part("imag" attribute) of the input function
	*	with a beta adapted to the simualtion betaOut. 
	*	The other attributes of jComponentIn are simply copied
	*
	* Description: 
	*   Reads the hybridization function from the jComponentIn array.
	* 	I name beta = jComponentIn["beta"]
	*	As we are in Fourrier space with Matsubara frequencies PI*(2*m+1)/beta, we need to adapt the read data to the beta in the simulation.
	*	This function computes function(iwn) from function(iwm) with wm = PI*(2*m+1)/beta, wn = PI*(2*n+1)/betaOut using linear interpolation.
	* 	The input function is (jComponentIn["real"] + i*jComponentIn["imag"])
	* 
	*/
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
		/** 
		* 
		* Function(std::string name, json const& jNumericalParams, json const& jEntry)
		* 
		* Parameters :	name : name of the Entry
		*				jNumericalParams : numerical parameters of the simulation
		*				jEntry : Entry
		*
		* Description: 
		*   Reads the matsubara entry from the json Object jEntry.
		*	Transforms the entry(Matsubara frequencies) to imaginary time with a Fourrier Transform
		* 	The little trick here is the following : 
		*		As the self-energy computed with this hybridation has high statisctical noise at high Matsubara frequency, we use a trick to compute the hybridation to correct this noise.
		*		We expand the hybridation function at the second order in i(omega_n) for omega->infinity and compute the Fourrier transform exactly
		*		We then compute numerically the remaining part of the fourrier transform thus reducing numerical noise.
		*		For more details, see appendix B.3 and B.1 of E. Gull's PHD Thesis (link: https://www.research-collection.ethz.ch/handle/20.500.11850/104013)
		* 	Stores the result for use in the Program
		*/
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
		/** 
		* get(double time)
		* 
		* Parameters :	time : time at wich we wànt the hybridation function
		*
		* Return Value : the  function at time `time` using linear interpolation on the discretized hybridation function we read from file
		*/
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
			jObject[el.key()] = read(jTemp[el.key()],betaOut);
		} 
	}
};

#endif