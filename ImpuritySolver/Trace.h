#ifndef __TRACE
#define __TRACE

#include <cmath>
#include <iostream>
#include <vector>
#include <climits>
#include <map>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <bitset>
#include <set>
#include "Utilities.h"
#include "MPIUtilities.h"
#include <boost/lexical_cast.hpp>


namespace Tr {
	/** 
	* 
	* struct Operator
	* 
	* Description: 
	*   This Class describes an operator (described by a time, a type). 
	*	The spin is not included here but directly in the bath (no diagonal term in the spin basis)
	*/
	struct Operator {
		Operator() : exp_(0) {};
		Operator(Operator const& other) : type_(other.type_), time_(other.time_), exp_(0) {};
		Operator(int type, double time) : type_(type), time_(time), exp_(0) {};
		Operator& operator=(Operator const& other) {
			type_ = other.type_;
			time_ = other.time_;
			
			return *this;
		};
		int type() const { return type_;};
		int* ptr() const { return 0;};
		double time() const { return time_;};
		/** 
		* 
		* std::complex<double> const* exp(unsigned int N, double beta)
		* 
		* Returns :	a vector of size N, his elements are exp(i*j*(time_*2*PI/beta)), j in [0,N-1]
		*
		* Description: 
		*   The type of the operator doesn't come in line
		*/
		std::complex<double> const* exp(unsigned int N, double beta) const {
			if(!exp_) {
				exp_ = new Ut::complex[N];
				
				Ut::complex const base = Ut::complex(std::cos(time_*2*M_PI/beta), std::sin(time_*2*M_PI/beta));
				Ut::complex entry = 1.;
				
				for(unsigned int n = 0; n < N; ++n) {
					exp_[n] = entry; entry *= base;
				}
			};
			return exp_;
		};
		
		~Operator() { 
			delete[] exp_;
		};
	private:
		int type_;
		double time_;
		mutable std::complex<double>* exp_;
	};

	
	bool operator<(Operator const& lhs, Operator const& rhs) {
		return lhs.time() < rhs.time();
	};
	
	struct Meas {
		Meas(json_spirit::mObject const& jNumericalParams) : 
		k(.0), N(.0), Sz(.0), D(.0), 
		Chi(.0,
			jNumericalParams.at("EObs").get_real() > .0 ? jNumericalParams.at("EObs").get_real()*jNumericalParams.at("beta").get_real()/(2.*M_PI) + 2 : 1
			) {
		};
		double k;
		double N;
		double Sz;
		double D;
		std::valarray<double> Chi;
	};
	/** 
	* 
	* struct Trace
	* 
	* Description: 
	*   This Class describes a site of the system. It stores all the informations about that site and allows all the updates.
	* 	The name Trace dosn't really mean anything in this program. Its name could be replaced by SiteHandler for example
	*/
	struct Trace {
		typedef std::set<Operator> Operators;
		/** 
		* 
		* Trace(json_spirit::mObject const& jNumericalParams, int site, Ut::Measurements& measurements,json_spirit::mObject const& jPreviousConfig)
		* 
		* Parameters :	jNumericalParams : storage of all the numerical parameters of the simulation
		*				site : number of the site this refers to
		*				measurements : variable used to store the measurements, used for output
		*				jPreviousConfig : previous configuration read from file
		*
		* Description: 
		*   If a previous configuration is available, insert all the indicated vertices. Else make it empty for the beginning
		*/
		Trace(json_spirit::mObject const& jNumericalParams, int site, Ut::Measurements& measurements,json_spirit::mObject const& jPreviousConfig) :
		beta_(jNumericalParams.at("beta").get_real()),
		U_(jNumericalParams.at("U").get_real()),
		mu_(jNumericalParams.at("mu").get_real()),
		shift_(std::max(0., std::max(mu_*beta_, 2.*mu_*beta_ - U_*beta_))),
		Tr00_0_(std::exp(-shift_)),
		Tr00_1_(std::exp(mu_*beta_ - shift_)),
		Tr00_2_(std::exp(2.*mu_*beta_ - U_*beta_ - shift_)),
		logTr00_(std::log(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_) + shift_),
		operators_(2),
		lenght_(2, .0),
		overlap_(.0),
		toChi_(0),
		acc_(jNumericalParams),
	    chiTemp_(acc_.Chi.size()) {
			operators_[0] = new Operators();
			operators_[1] = new Operators();
			
			if(jPreviousConfig.size()) {
				json_spirit::mObject const& jPreviousConfigSite = jPreviousConfig.at("Site " + boost::lexical_cast<std::string>(site)).get_obj();
				for(int spin = 0; spin < 2; ++spin) {
					json_spirit::mArray const&  jPreviousConfigSiteSpin = jPreviousConfigSite.at("Spin " + boost::lexical_cast<std::string>(spin)).get_array();
					std::size_t size = jPreviousConfigSiteSpin.size();
					for(std::size_t i = 0; i < size; ++i) {
						int type = jPreviousConfigSiteSpin[i].get_obj().at("type").get_int();
						double time = jPreviousConfigSiteSpin[i].get_obj().at("time").get_real();
						operators_[spin]->insert(Operator(type, time));
					}
				}
				set();
			}
		};
		
		Operator const& op() { return op_;};
		Operator const& opDagg() { return opDagg_;};
		
		Operators const& operators(int spin) { return *operators_[spin];};
		
		int sign() { 
			int sign = 1;
			if(operators(0).size()%4 && !operators(0).begin()->type()) sign *= -1;
			if(operators(1).size()%4 && !operators(1).begin()->type()) sign *= -1;
			return sign;
		};
		/** 
		* 
		* double insert(int spin, Ut::UniformRngType& urng)
		* 
		* Parameters :	spin : spin of the new vertex
		*				urng : probabiliy distribution
		*
		* Return Value : the Metropolis Hastings acceptation rate for the site part 
		*
		* Description: 
		*   In this function we want to add a new vertex. 
		* 	To do so, we choose what type of vertex and where to insert it.
		*	We then insert it and compute the part of the acceptation rate relative to this site only and return it.
		*/
		double insert(int spin, Ut::UniformRngType& urng) {
			Operators& ops = *operators_[spin]; Operators& opsOther = *operators_[1 - spin]; 
			
			std::pair<Operators::iterator, bool> insertInfo;
			//Insert one new vertex in ops
			do insertInfo = ops.insert(Operator(0, beta_*urng())); while(!insertInfo.second);
			//get the inserted vertex
			op_ = *insertInfo.first; 
			
			double l;
			Operators::iterator itTemp = insertInfo.first;
			Operator opNext = itTemp == ops.begin() ? *ops.rbegin() : *--itTemp;
			if(!opNext.type()) {
				l = control(op_.time() - opNext.time());
				do insertInfo = ops.insert(Operator(1, control(op_.time() - urng()*l))); while(!insertInfo.second);
				opDagg_ = *insertInfo.first;
				
				lenghtDiff_ = control(op_.time() - opDagg_.time());
				overlapDiff_ = opDagg_ < op_ ? otherLenght(opsOther, opDagg_, op_) : otherLenghtWO(opsOther, opDagg_, op_);
			} else {
				opNext = *(++insertInfo.first == ops.end() ? ops.begin() : insertInfo.first);
				
				l = control(opNext.time() - op_.time());
				do insertInfo = ops.insert(Operator(1, control(opNext.time() - urng()*l))); while(!insertInfo.second);
				opDagg_ = *insertInfo.first;
				
				lenghtDiff_ = -control(opDagg_.time() - op_.time());
				overlapDiff_ = -(op_ < opDagg_ ? otherLenght(opsOther, op_, opDagg_) : otherLenghtWO(opsOther, op_, opDagg_));
			}
			
			double const qRatio = std::log(beta_*l/(ops.size() > 2 ? ops.size() : 1));

			//On retourne la probabilité d'acceptation suivant les différents cas particuliers
			if(ops.size() > 2 && opsOther.size()) return lenghtDiff_*mu_ - overlapDiff_*U_ + qRatio;
			if(ops.size() > 2) return logTr0(lenght_[spin] + lenghtDiff_) - logTr0(lenght_[spin]) + qRatio;
		    if(opsOther.size()) return (lenght_[1 - spin] + lenghtDiff_)*mu_ - overlapDiff_*U_ - logTr0(lenght_[1 - spin]) + qRatio; 
			return logTr0(lenghtDiff_) - logTr00_ + qRatio;
		};
		/** 
		* 
		* int acceptInsert(int spin)
		* 
		* Parameters :	spin : spin of the new vertex
		*
		* Return Value : the sign change according to the ordering of the operators 
		*
		* Description: 
		*	It updates the overlap and occupation observables for the site.
		*/
		int acceptInsert(int spin) {
			lenght_[spin] += lenghtDiff_;
			overlap_ += overlapDiff_;
			
			delete[] toChi_; toChi_ = 0;
			
			return op_ < opDagg_ ? -1 : 1;
		};
		/** 
		* 
		* void rejectInsert(int spin)
		* 
		* Parameters :	spin : spin of the new vertex
		*
		* Description: 
		*	It removes the operators we inserted in the insert function, as the update was rejected.
		*/
		
		void rejectInsert(int spin) {
			operators_[spin]->erase(op_);
			operators_[spin]->erase(opDagg_);
		};
		/** 
		* 
		* double erase(int spin, Ut::UniformRngType& urng)
		* 
		* Parameters :	spin : spin of the vertex to remove
		*				urng : probabiliy distribution
		*
		* Return Value : the Metropolis Hastings acceptation rate for the site part 
		*
		* Description: 
		*   In this function we want to remove a vertex. 
		* 	To do so, we choose what vertex we want to remove (accoring to urng).
		*	We then compute the part of the acceptation rate relative to this site only and return it.
		*/
		double erase(int spin, Ut::UniformRngType& urng) {
			Operators const& ops = operators(spin); Operators const& opsOther = operators(1 - spin); 
			
			if(ops.size() != 2) {
				Operators::iterator itLow = ops.begin();
				std::advance(itLow, static_cast<int>(ops.size()*urng())); 
				Operators::iterator itUp = itLow; 
				if(++itUp == ops.end()) itUp = ops.begin(); 
				
				double l;
				if(itLow->type()) {
					op_ = *itUp; opDagg_ = *itLow;
					lenghtDiff_ = -control(itUp->time() - itLow->time());
					overlapDiff_ = -(*itLow < *itUp ? otherLenght(opsOther, *itLow, *itUp) : otherLenghtWO(opsOther, *itLow, *itUp));
					if(itLow == ops.begin()) itLow = ops.end();
					l = control(itUp->time() - (--itLow)->time());
				} else { 
					op_ = *itLow; opDagg_ = *itUp;
					lenghtDiff_ = control(itUp->time() - itLow->time());
					overlapDiff_ = *itLow < *itUp ? otherLenght(opsOther, *itLow, *itUp) : otherLenghtWO(opsOther, *itLow, *itUp);
					if(++itUp == ops.end()) itUp = ops.begin();
					l = control(itUp->time() - itLow->time());
				}
				
				double const qRatio = std::log(ops.size()/(beta_*l));
				
				if(opsOther.size()) return lenghtDiff_*mu_ - overlapDiff_*U_ + qRatio;
				return logTr0(lenght_[spin] + lenghtDiff_) - logTr0(lenght_[spin]) + qRatio;
			} else {
				op_ = *ops.begin(); opDagg_ = *ops.rbegin(); if(op_.type()) std::swap(op_, opDagg_);
				lenghtDiff_ = -lenght_[spin]; overlapDiff_ = -overlap_; 
				
				double const qRatio = -std::log(beta_*beta_);
				
				if(opsOther.size()) return logTr0(lenght_[1 - spin]) - ((lenght_[spin] + lenght_[1 - spin])*mu_ - U_*overlap_) + qRatio; 
				return logTr00_ - logTr0(lenght_[spin]) + qRatio;		
			} 
		};
		/** 
		* 
		* int acceptErase(int spin)
		* 
		* Parameters :	spin : spin of the vertex to remove
		*
		* Return Value : the sign change according to the ordering of the operators 
		*
		* Description: 
		*	It updates the overlap and occupation observables for the site.
		*/
		int acceptErase(int spin) {
			lenght_[spin] += lenghtDiff_; 
			overlap_ += overlapDiff_;
			
			operators_[spin]->erase(op_);
			operators_[spin]->erase(opDagg_);
			
			delete[] toChi_; toChi_ = 0;
			
			return op_ < opDagg_ ? -1 : 1;
		};
		
		void rejectErase(int spin) {			
		};
		
		/** 
		* 
		* void flip()
		*
		* Description: 
		*	It flips the operators for the two spins and also the lenght_ which is directly a function of the operators
		*/
		void flip() {
			std::swap(lenght_[0], lenght_[1]);
			std::swap(operators_[0], operators_[1]);
		};

		/** 
		* 
		* measure(int sign) 
		* 
		* Parameters :	sign : the sign of the system for the obersvables to be printed correctly
		*
		* Description: 
		*   Computes all the observables for the site (Occupation, Double Occupation, Spin, Chi) 
		* 
		*/
		
		std::valarray<Ut::complex>& measure(int sign) {
			chiTemp_ = 0; //We reset the temporary Chiij to make sure we don't make mistakes later
			acc_.k += sign*((operators(0).size() + operators(1).size())/2.);
			
			//Here we compute the observables depending on the number of operators we have in our timeline
			if(operators(0).size() && operators(1).size()) {
				acc_.N += sign*.5*(lenght_[0] + lenght_[1])/beta_;
				acc_.D += sign*overlap_/beta_;
				acc_.Sz += sign*.5*(lenght_[0] - lenght_[1])/beta_;
				acc_.Chi[0] += sign*.25*(lenght_[0] - lenght_[1])*(lenght_[0] - lenght_[1])/beta_;

				//We want the two-site dependant spin susceptibility, we get the Chiij for each site and multiply them later.
				chiTemp_[0] = .5*(lenght_[0] - lenght_[1]);
			} else if(operators(0).size()) {
				double arg = mu_*beta_ - U_*lenght_[0];
				double exp = std::exp(-std::abs(arg));
				double fact = (arg < .0 ? exp/(1. + exp) : 1./(1. + exp));
				
				acc_.N += sign*(beta_*fact + lenght_[0])/(2.*beta_);
				acc_.D += sign*lenght_[0]*fact/beta_;
				acc_.Sz += sign*(-beta_*fact + lenght_[0])/(2.*beta_);
				acc_.Chi[0] += sign*.25*((beta_ - 2.*lenght_[0])*fact + lenght_[0]*lenght_[0]/beta_);
				//this next expression is only good when for different sites. This is taken into account in the Markovchain file
				chiTemp_[0] = .5*(-beta_*fact + lenght_[0]);
			} else if(operators(1).size()) {
				double arg = mu_*beta_ - U_*lenght_[1];
				double exp = std::exp(-std::abs(arg));
				double fact = (arg < .0 ? exp/(1. + exp) : 1./(1. + exp));
				
				acc_.N += sign*(beta_*fact + lenght_[1])/(2.*beta_);
				acc_.D += sign*lenght_[1]*fact/beta_;
				acc_.Sz += sign*(beta_*fact - lenght_[1])/(2.*beta_);
				acc_.Chi[0] += sign*.25*((beta_ - 2.*lenght_[1])*fact + lenght_[1]*lenght_[1]/beta_);
				//This expression is only good when for different sites. This is taken into account in the Markovchain file
				chiTemp_[0] =  .5*(beta_*fact - lenght_[1]);
			} else {
				acc_.N += sign*(Tr00_1_ + Tr00_2_)/(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_);
				acc_.D += sign*Tr00_2_/(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_);
				acc_.Sz = .0;
				acc_.Chi[0] += sign*beta_/2.*Tr00_1_/(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_);
				//This expression is only good when for different sites. This is taken into account in the Markovchain file
				chiTemp_[0] += .0;
			}
			

			if(acc_.Chi.size() > 1) {
				if(!toChi_) {     //So wies aussieht haben wir hier glueck: keine fall unterschiedung fŸr 0 expansions ordnung nštig fuer finites omega	
					//This function is doing the same job as if it were only going through all the operators and calculating chi					
					toChi_ = new Ut::complex[acc_.Chi.size()];
					
					Operators::const_iterator it0 = operators(0).begin();
					Operators::const_iterator it1 = operators(1).begin();
					
					while(it0 != operators(0).end() || it1 != operators(1).end()) {
						double const time1 = it1 != operators(1).end() ? it1->time() : std::numeric_limits<double>::max();
						
						while(it0 != operators(0).end() && it0->time() < time1) {	
							int const one = 1; int const N = acc_.Chi.size(); Ut::complex const DSz = .5*(1 - 2*it0->type());
							zaxpy_(&N, &DSz, it0->exp(N, beta_), &one, toChi_, &one);
							it0++;					 
						}
						
						double const time0 = it0 != operators(0).end() ? it0->time() : std::numeric_limits<double>::max();
						
						while(it1 != operators(1).end() && !(time0 < it1->time())) {
							int const one = 1; int const N = acc_.Chi.size(); Ut::complex const DSz = -.5*(1 - 2*it1->type());
							zaxpy_(&N, &DSz, it1->exp(N, beta_), &one, toChi_, &one);
							it1++;
						}
					}
				}

				for(unsigned int n = 1; n < acc_.Chi.size(); ++n){
					chiTemp_[n] = toChi_[n]; 
					acc_.Chi[n] += sign*(toChi_[n].real()*toChi_[n].real() + toChi_[n].imag()*toChi_[n].imag());				
				}
			};
			return chiTemp_;
		};
		/** 
		* void store(Ut::Measurements& measurements, int site, Meas& meas, int NMeas)
		* 
		* Parameters :	measurements : variable used to store the measurements, used for output
		*				site : number of the site (because we don't store it in the Object)
		*				meas : measurement object of the entire simulation
		*				NMeas : Number of measurements done since the last time we stored some measurements
		*
		* Description: 
		* 	Computes the mean over all measurements done since the last time we stored them.
		*   Stores the measurements for this site in the handler.
		*	Adds the site measurements to the measurements of the entire simulation
		*	Resets the measurements for the next round of measurements
		*/
		void store(Ut::Measurements& measurements, int site, Meas& meas, int NMeas) {
			acc_.k /= NMeas; 
			acc_.N /= NMeas; 
			acc_.Sz /= NMeas; 
			acc_.D /= NMeas; 
			acc_.Chi /= NMeas;
			
			std::string s = boost::lexical_cast<std::string>(site);
			measurements["k_" + s] << acc_.k;
			measurements["N_" + s] << acc_.N;
			measurements["D_" + s] << acc_.D;
			measurements["Sz_" + s] << acc_.Sz;
			measurements["Chi0_" + s] << acc_.Chi[0];
			
			if(acc_.Chi.size() > 1) {
				for(unsigned int n = 1; n < acc_.Chi.size(); ++n) 
					acc_.Chi[n] *= beta_/((2*n*M_PI)*(2*n*M_PI)); 
				
				measurements["Chi_" + s] << acc_.Chi;
			};
			
			meas.k += acc_.k; acc_.k = .0;
			meas.N += acc_.N; acc_.N = .0;
			meas.Sz += acc_.Sz; acc_.Sz = .0;
			meas.D += acc_.D; acc_.D = .0;
			meas.Chi += acc_.Chi; acc_.Chi = .0;
		};	
		/** 
		* std::valarray<double>& getChi()
		* 
		* No Parameters
		*
		* Description: 
		* 	Returns the Chi accumulator. This is needed in the MarkovChain.h file in order to compute correctly the Spin Susceptibility because of the terms
		*	with no electrons on the segments (as seen in the measure function). See Readme for more explanation
		*/
		std::valarray<double>& getChi(){
			return acc_.Chi;
		}	
		/** 
		* void saveConfig(json_spirit::mObject& jConfig,int site)
		* 
		* Parameters :	jConfig : Object in which the configuration should be saved
		*				site : number of the site (because we don't store it in the Object)
		*
		* Description: 
		* 	Saves the current simulation state (the segement operators) to the jConfig object.
		*/
		void saveConfig(json_spirit::mObject& jConfig,int site) {
			json_spirit::mObject jConfigSite;
			for(int spin = 0; spin < 2; ++spin) {
				json_spirit::mArray jConfigSiteSpin;
				for(Operators::const_iterator it = operators(spin).begin();  it != operators(spin).end(); ++it) {
					json_spirit::mObject jElement;
					jElement["type"] = it->type();
					jElement["time"] = it->time();
					jConfigSiteSpin.push_back(jElement);
				}
				jConfigSite["Spin " + boost::lexical_cast<std::string>(spin)] = jConfigSiteSpin;
			}
			jConfig["Site "+ boost::lexical_cast<std::string>(site)] = jConfigSite;
		};
		
		~Trace() {
			delete[] toChi_;
			
			delete operators_[0];
			delete operators_[1];
 		};
	private:
		double const beta_;
		double const U_;
		double const mu_;
		
		double const shift_;
		double const Tr00_0_;
		double const Tr00_1_;
		double const Tr00_2_;
		double const logTr00_;
		
		std::vector<Operators*> operators_;
		
		std::vector<double> lenght_;
		double overlap_;
		
		Ut::complex* toChi_;
		Meas acc_;		
		std::valarray<Ut::complex> chiTemp_;
		
		double lenghtDiff_;
		double overlapDiff_;
		
		Operator op_;
		Operator opDagg_;
		
		inline double control(double tau) {
			return tau <= .0 ? tau + beta_ : tau;
		};
		
		/** 
		* double otherLenght(Operators const& opsOther, Operator const& opLow, Operator const& opUp)
		* 
		* Parameters :	opsOther : all the operators for the opposite spin
		*				opLow : operator with the lowest time
		*			 	opUp : operator with the highest time
		*
		* Return Value : The overlap time between this vertex and all vertex of the opsOther table
		*
		* Description: 
		* 	We assume 0 < opLow.time() < opUp.time() < beta
		* 	Iterates over all operators in OpsOther to find the Double occupation associated with this pair of operators
		*/
		double otherLenght(Operators const& opsOther, Operator const& opLow, Operator const& opUp) { 
			if(opsOther.size()) {
				Operators::const_iterator it = opsOther.upper_bound(opLow);
				if(it != opsOther.end()) {
					int temp = 1 - 2*it->type();
					
					double over = -.5*(1 + temp)*opLow.time();
					while(it != opsOther.end() && *it < opUp) {
						over += temp*(it++)->time();
						temp *= -1;
					}
					over += .5*(1 + temp)*opUp.time();
					
					return over;
				} else 
					return (--it)->type()*(opUp.time() - opLow.time());
			} else 
				return .0;
		};
		/** 
		* double otherLenght(Operators const& opsOther, Operator const& opLow, Operator const& opUp)
		* 
		* Parameters :	opsOther : all the operators for the opposite spin
		*				opLow : operator with the highest time
		*			 	opUp : operator with the lowest time
		*
		* Return Value : The overlap time between this vertex and all vertex of the opsOther table
		*
		* Description:
		* 	We assume 0 < opUp.time() < opLow.time() < beta
		*	This has the same effect as otherLenght above, except the operators are inverted.
		* 	So we compute the overlap from opLow.time() to beta an then from 0 to opUp.time()
		*/
		double otherLenghtWO(Operators const& opsOther, Operator const& opLow, Operator const& opUp) {
			if(opsOther.size()) {
				Operators::const_iterator it = opsOther.upper_bound(opLow);
				int temp = (1 - 2*(it != opsOther.end() ? it->type() : opsOther.begin()->type()));
				
				double over = -.5*(1 + temp)*opLow.time();
				while(it != opsOther.end()) {
					over += temp*(it++)->time();
					temp *= -1;
				}       
				it = opsOther.begin();
				while(it != opsOther.end() && *it < opUp) {
					over += temp*((it++)->time() + beta_);
					temp *= -1;
				}
				over += .5*(1 + temp)*(opUp.time() + beta_);
				
				return over;
			} else 
				return .0;
		};
		
		double logTr0(double l) {
			double const arg = beta_*mu_ - U_*l;
			double const exp = std::exp(-std::abs(arg));
			return l*mu_ + std::log(1. + exp) + (arg > .0 ? arg : .0);
		};
		
		/** 
		* void set()
		*
		* Description:
		* 	This function computes the occupation length and the overlap directly from the operators. (instead of the usual step by step approach)
		* 	This can be used at loading time to start the simulation or whenever is needed to recompute from scratch the overlap and the occupation
		*/
		void set() {
			std::set<std::pair<Operator, int> > temp;
			std::vector<int> isSeg(2);
			
			for(int spin = 0; spin < 2; ++spin)
			{ 
				if(operators(spin).size()) {
					isSeg[spin] = !operators(spin).begin()->type();
					
					lenght_[spin] = !operators(spin).begin()->type() ? beta_ : .0;
					for(Operators::iterator it = operators(spin).begin(); it != operators(spin).end(); ++it) {
						lenght_[spin] += it->type() ? -it->time() : it->time();
						temp.insert(std::make_pair(*it, spin));
					}
				} else 
					lenght_[spin] = .0;
			}
			
			if(operators(0).size() && operators(1).size()) {
				overlap_ = isSeg[0] && isSeg[1] ? beta_ : .0;
				for(std::set<std::pair<Operator, int> >::iterator it = temp.begin(); it != temp.end(); ++it) {
					if(isSeg[1 - it->second])
					{
						overlap_ += it->first.type() ? -it->first.time() : it->first.time();
					} 
					isSeg[it->second] = it->first.type();
				}
			} else 
				overlap_ = .0; 
			
		};
	};
};

#endif






















