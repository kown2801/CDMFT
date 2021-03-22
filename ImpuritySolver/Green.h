#ifndef __GREEN
#define __GREEN

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <sstream>
#include "nlohmann_json.hpp"
#include "Utilities.h"

namespace Green {
	struct Meas {
		Meas(std::string name, json const& jNumericalParams, Ut::Measurements& measurements) :
		name_(name),
		beta_(jNumericalParams["beta"]),
		nMatG_(beta_*jNumericalParams["EGreen"].get<double>()/(2*M_PI) + 1), 
		nItG_(4*(2*nMatG_ + 1)), 
		DeltaInv_(nItG_/beta_),
		green_(new double[4*nItG_]) {
			std::memset(green_, 0, 4*nItG_*sizeof(double));	
		};
		
		void add(double time, double value) {
			int index = static_cast<int>(DeltaInv_*time);
			double Dtime = time - static_cast<double>(index + .5)/DeltaInv_;
			
			double* green = green_ + 4*index;
			*green++ += value; 
			value *= Dtime;
			*green++ += value;
			value *= Dtime;
			*green++ += value;
			value *= Dtime;
			*green += value;
		};
		
		void measure(Ut::Measurements& measurements, int NAlpsMeas) {
			std::valarray<double> greenReal(nMatG_);
			std::valarray<double> greenImag(nMatG_);
			double Dtau = beta_/static_cast<double>(nItG_);
			
			for(int m = 0; m < nMatG_; ++m) {
				double omega = M_PI*static_cast<double>(2*m + 1)/beta_;
				double lambda = -2.*std::sin(omega*Dtau/2.)/((Dtau*omega*(1. - omega*omega*Dtau*Dtau/24.))*beta_*NAlpsMeas);  //missing -1/beta factor
				
				Ut::complex iomega(.0, omega);
				Ut::complex exp(std::exp(iomega*Dtau/2.));
				Ut::complex fact(std::exp(iomega*Dtau));  
				
				Ut::complex tempMAT(.0);
				double const* itGreen = green_;
				for(int i = 0; i < nItG_; ++i) {
					Ut::complex coeff = lambda*exp;
					tempMAT += coeff**itGreen++;
					coeff *= iomega;
					tempMAT += coeff**itGreen++;
					coeff *= iomega/2.;
					tempMAT += coeff**itGreen++;
					coeff *= iomega/3.;
					tempMAT += coeff**itGreen++;
					
					exp *= fact;
				}
				
				greenReal[m] = tempMAT.real();
				greenImag[m] = tempMAT.imag();
			}
			
			measurements["GreenR_" + name_] << greenReal;
			measurements["GreenI_" + name_] << greenImag; 
			
			std::memset(green_, 0, 4*nItG_*sizeof(double));
		};
		~Meas() { delete[] green_;};
	private:
		std::string const name_;
		
		double const beta_;
		
		int const nMatG_;
		int const nItG_;
		double const DeltaInv_;
		
		double* const green_;
	};
	
};
#endif