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
#include <set>
#include "Utilities.h"
#include "MPIUtilities.h"
#include <boost/lexical_cast.hpp>


/////////////////////////////////////////////////////////////////
//Scheiss Observables !!!!!!!
// statistik für k auch individuell für \up und \down möglich
/////////////////////////////////////////////////////////////////

namespace Tr {
	//Da beim erase die operatoren nicht veraendert werden gehen keine berechneten exponentiale verlohren
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
	
	struct Trace {
		typedef std::set<Operator> Operators;
		
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
		acc_(jNumericalParams) {
			operators_[0] = new Operators();
			operators_[1] = new Operators();
			
			if(jPreviousConfig.size()) {
				json_spirit::mObject jPreviousConfigSite = jPreviousConfig.at("Site " + boost::lexical_cast<std::string>(site)).get_obj();
				for(int spin = 0; spin < 2; ++spin) {
					json_spirit::mArray jPreviousConfigSiteSpin = jPreviousConfigSite.at("Spin " + boost::lexical_cast<std::string>(spin)).get_array();
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
			if(operators(1).size()%4 &&  operators(1).begin()->type()) sign *= -1;
			return sign;
		};
		
		double insert(int spin, Ut::UniformRngType& urng) {
			Operators& ops = *operators_[spin]; Operators& opsOther = *operators_[1 - spin]; 
			
			std::pair<Operators::iterator, bool> insertInfo;
			do insertInfo = ops.insert(Operator(0, beta_*urng())); while(!insertInfo.second);
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
			
			if(ops.size() > 2 && opsOther.size()) return lenghtDiff_*mu_ - overlapDiff_*U_ + qRatio;
			if(ops.size() > 2) return logTr0(lenght_[spin] + lenghtDiff_) - logTr0(lenght_[spin]) + qRatio;
		    if(opsOther.size()) return (lenght_[1 - spin] + lenghtDiff_)*mu_ - overlapDiff_*U_ - logTr0(lenght_[1 - spin]) + qRatio; 
			return logTr0(lenghtDiff_) - logTr00_ + qRatio;
		};
		
		int acceptInsert(int spin) {
			lenght_[spin] += lenghtDiff_;
			overlap_ += overlapDiff_;
			
			delete[] toChi_; toChi_ = 0;
			
			return op_ < opDagg_ ? 2*spin - 1 : 1 - 2*spin;
		};
		
		void rejectInsert(int spin) {
			operators_[spin]->erase(op_);
			operators_[spin]->erase(opDagg_);
		};
		
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
	
		int acceptErase(int spin) {
			lenght_[spin] += lenghtDiff_; 
			overlap_ += overlapDiff_;
			
			operators_[spin]->erase(op_);
			operators_[spin]->erase(opDagg_);
			
			delete[] toChi_; toChi_ = 0;
			
			return op_ < opDagg_ ? 2*spin - 1 : 1 - 2*spin;
		};
		
		void rejectErase(int spin) {			
		};
		
		void flip() {
			std::swap(lenght_[0], lenght_[1]);
			std::swap(operators_[0], operators_[1]);
		};
		
		void measure(int sign) {
			acc_.k += sign*((operators(0).size() + operators(1).size())/2.);
			
			if(operators(0).size() && operators(1).size()) {
				acc_.N += sign*.5*(lenght_[0] + lenght_[1])/beta_;
				acc_.D += sign*overlap_/beta_;
				acc_.Sz += sign*.5*(lenght_[0] - lenght_[1])/beta_;
				acc_.Chi[0] += sign*.25*(lenght_[0] - lenght_[1])*(lenght_[0] - lenght_[1])/beta_;
			} else if(operators(0).size()) {
				double arg = mu_*beta_ - U_*lenght_[0];
				double exp = std::exp(-std::abs(arg));
				double fact = (arg < .0 ? exp/(1. + exp) : 1./(1. + exp));
				
				acc_.N += sign*(beta_*fact + lenght_[0])/(2.*beta_);
				acc_.D += sign*lenght_[0]*fact/beta_;
				acc_.Sz += sign*(-beta_*fact + lenght_[0])/(2.*beta_);
				acc_.Chi[0] += sign*.25*((beta_ - 2.*lenght_[0])*fact + lenght_[0]*lenght_[0]/beta_);
			} else if(operators(1).size()) {
				double arg = mu_*beta_ - U_*lenght_[1];
				double exp = std::exp(-std::abs(arg));
				double fact = (arg < .0 ? exp/(1. + exp) : 1./(1. + exp));
				
				acc_.N += sign*(beta_*fact + lenght_[1])/(2.*beta_);
				acc_.D += sign*lenght_[1]*fact/beta_;
				acc_.Sz += sign*(beta_*fact - lenght_[1])/(2.*beta_);
				acc_.Chi[0] += sign*.25*((beta_ - 2.*lenght_[1])*fact + lenght_[1]*lenght_[1]/beta_);
			} else {
				acc_.N += sign*(Tr00_1_ + Tr00_2_)/(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_);
				acc_.D += sign*Tr00_2_/(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_);
				acc_.Sz = .0;
				acc_.Chi[0] += sign*beta_/2.*Tr00_1_/(Tr00_0_ + 2.*Tr00_1_ + Tr00_2_);
			}
			
			if(acc_.Chi.size() > 1) {
				if(!toChi_) {     //So wies aussieht haben wir hier glueck: keine fall unterschiedung für 0 expansions ordnung nötig fuer finites omega					
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

				for(unsigned int n = 1; n < acc_.Chi.size(); ++n) 
					acc_.Chi[n] += sign*(toChi_[n].real()*toChi_[n].real() + toChi_[n].imag()*toChi_[n].imag());				
			};
		};
		
		void measure(Ut::Measurements& measurements, int site, Meas& meas, int NAlpsMeas) {
			acc_.k /= NAlpsMeas; 
			acc_.N /= NAlpsMeas; 
			acc_.Sz /= NAlpsMeas; 
			acc_.D /= NAlpsMeas; 
			acc_.Chi /= NAlpsMeas;
			
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
			}
			
			meas.k += acc_.k; acc_.k = .0;
			meas.N += acc_.N; acc_.N = .0;
			meas.Sz += acc_.Sz; acc_.Sz = .0;
			meas.D += acc_.D; acc_.D = .0;
			meas.Chi += acc_.Chi; acc_.Chi = .0;
		}
		
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
		}
		
		/*
		void test() {
			std::vector<double> lenghtTemp(lenght_);
			double overlapTemp = overlap_;
			set();
			
			if(std::abs(lenght_[0] - lenghtTemp[0]) > 1.e-12 || std::abs(lenght_[0] - lenghtTemp[0]) > 1.e-12) 
				std::cout << "problem mit lenght" << std::endl;
			
			if(std::abs(overlap_ - overlapTemp) > 1.e-12)
				std::cout << "problem mit overlap" << std::endl;
		};
		*/

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
		
		double lenghtDiff_;
		double overlapDiff_;
		
		Operator op_;
		Operator opDagg_;
		
		inline double control(double tau) {
			return tau <= .0 ? tau + beta_ : tau;
		};
		
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
		
		void set() {
			std::set<std::pair<Operator, int> > temp;
			std::vector<int> isSeg(2);
			
			for(int spin = 0; spin < 2; ++spin) 
				if(operators(spin).size()) {
					isSeg[spin] = !operators(spin).begin()->type();
					
					lenght_[spin] = !operators(spin).begin()->type() ? beta_ : .0;
					for(Operators::iterator it = operators(spin).begin(); it != operators(spin).end(); ++it) {
						lenght_[spin] += it->type() ? -it->time() : it->time();
						temp.insert(std::make_pair(*it, spin));
					}
				} else 
					lenght_[spin] = .0;
			
			if(operators(0).size() && operators(1).size()) {
				overlap_ = isSeg[0] && isSeg[1] ? beta_ : .0;
				for(std::set<std::pair<Operator, int> >::iterator it = temp.begin(); it != temp.end(); ++it) {
					if(isSeg[1 - it->second]) overlap_ += it->first.type() ? -it->first.time() : it->first.time();
					isSeg[it->second] = it->first.type();
				}
			} else 
				overlap_ = .0; 
			
		};
	};
};

#endif






















