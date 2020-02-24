#ifndef __BATH
#define __BATH

#include <cmath>
#include <iostream>
#include <vector>
#include <climits>
#include <map>
#include <algorithm>
#include <fstream>
#include <cassert>
#include "Utilities.h"

//!! Achtung mit der Zeitordnung !!

namespace Ba {
	struct Matrix {
		Matrix(int I, int J) : I_(I), J_(J), data_(new double[I*J]) {};
		double& at(int i, int j) { return data_[i + I_*j];}; 
		double const& at(int i, int j) const { return data_[i + I_*j];};
		double* data() { return data_;};
		double const* data() const { return data_;};
		double* data(int i, int j) { return data_ + i + j*I_;};
		double const* data(int i, int j) const { return data_ + i + j*I_;};
		~Matrix() { delete[] data_;}
	private:
		int const I_;
		int const J_;
		double* const data_;
	};
	
	struct Operator {
		Operator() {};
		Operator(int site, int spin, double time, int* ptr) : site_(site), spin_(spin), time_(time), ptr_(ptr) {};
		int site() const { return site_;};
		int spin() const { return spin_;};
		double time() const { return time_;};
		int n() const { return *ptr_;};
	private:
		int site_;
		int spin_;
		double time_;
		int* ptr_;
	};
	
	bool operator==(Operator const& lhs, Operator const& rhs) {
		return lhs.site() == rhs.site() && lhs.spin() == rhs.spin() && lhs.time() == rhs.time();
	}
	
	struct Bath {
		struct GreenIterator {
			GreenIterator(Operator const* opBegin, Operator const* opEnd, Operator const* opDagg, double const* value) : opBegin_(opBegin), opEnd_(opEnd), opR_(opBegin), opL_(opDagg), value_(value) {};
			Operator const& opR() const { return *opR_;};
			Operator const& opL() const { return *opL_;};			
			double value() const { return *value_;};
			
			GreenIterator& operator++() { ++value_; if(++opR_ == opEnd_) { opR_ = opBegin_; ++opL_;}; return *this;}; 
			bool operator!=(GreenIterator const& other) const { return value_ != other.value_;};
		private:
			Operator const* const opBegin_; Operator const* const opEnd_;
			Operator const* opR_; Operator const* opL_; double const* value_;
		};
		
		Bath() : swapSign_(1), B_(0), det_(1.) {};
		
		void add(int site, int spin, double time, int* ptr) { spin ? opsL_.push_back(Operator(site, spin, time, ptr)) : opsR_.push_back(Operator(site, spin, time, ptr));};
		void addDagg(int site, int spin, double time, int* ptr) { spin ? opsR_.push_back(Operator(site, spin, time, ptr)) : opsL_.push_back(Operator(site, spin, time, ptr));};
		
		GreenIterator begin() const { return B_ ? GreenIterator(opsR_.data(), opsR_.data() + opsR_.size(), opsL_.data(), B_->data()) : GreenIterator(0, 0, 0, 0);};
		GreenIterator end() const { return B_ ? GreenIterator(0, 0, 0, B_->data() + opsR_.size()*opsL_.size()) : GreenIterator(0, 0, 0, 0);};
		
		//template<typename T, typename L> Bath(int spin, std::vector<Tr*>::const_iterator begin, std::vector<Tr*>::const_iterator end, L const& link) {
		//};
		
		template<typename L> double insert(int site, int spin, double time, int* ptr, double timeDagg, int* ptrDagg, L const& link) {
			opR_ = Operator(site, spin, time, ptr); 
			opL_ = Operator(site, spin, timeDagg, ptrDagg);
			if(spin) std::swap(opR_, opL_);
			
			int const N = opsL_.size();
			val_ = link(opL_, opR_);
			
			if(N) {  
				Bv_.resize(N); vec_.resize(N);
				for(int n = 0; n < N; ++n) vec_[n] = link(opsL_[n], opR_);
				
				char const no = 'n';
				int const inc = 1;
				double const zero = .0;
				double const one = 1.;
				dgemv_(&no, &N, &N, &one, B_->data(), &N, vec_.data(), &inc, &zero, Bv_.data(), &inc);
				
				for(int n = 0; n < N; ++n) vec_[n] = link(opL_, opsR_[n]);			
				val_ -= ddot_(&N, vec_.data(), &inc, Bv_.data(), &inc);	
			}
			
			return std::log(std::abs(val_));
		}; 
		
		int acceptInsert() {
			int const N = opsL_.size();
			int const newN = N + 1;
			
			opsR_.push_back(opR_); opsL_.push_back(opL_); 
			
			Matrix* temp = new Matrix(newN, newN);
			
			double fact = 1./val_;
			temp->at(N, N) = fact;
			
			if(N) {                                    
				std::vector<double> hBTilde(N);
				
				char const yes = 't';
				int const inc = 1;
				double const zero = .0;
				double const one = 1.;
				dgemv_(&yes, &N, &N, &fact, B_->data(), &N, vec_.data(), &inc, &zero, hBTilde.data(), &inc);
				dger_(&N, &N, &one, Bv_.data(), &inc, hBTilde.data(), &inc, B_->data(), &N);
				
				for(int n = 0; n < N; ++n) 
					dcopy_(&N, B_->data(0, n), &inc, temp->data(0, n), &inc);
				
				fact = -1./val_;
				dscal_(&N, &fact, Bv_.data(), &inc); 
				dcopy_(&N, Bv_.data(), &inc, temp->data(0, N), &inc);
				
				double const minus = -1.;
				dscal_(&N, &minus, hBTilde.data(), &inc); 
				dcopy_(&N, hBTilde.data(), &inc, temp->data(N, 0), &newN); 
			}
			
			delete B_; B_ = temp;
			
			det_ *= val_;
			
			return val_ > .0 ? 1 : -1;
		};
		
		double erase(int site, int spin, double time, double timeDagg) {
			Operator dummyR(site, spin, time, 0);
			Operator dummyL(site, spin, timeDagg, 0);
			if(spin) std::swap(dummyR, dummyL);
			   
			posR_ = 0; while(static_cast<std::size_t>(posR_) != opsR_.size() && !(opsR_[posR_] == dummyR)) ++posR_;
			posL_ = 0; while(static_cast<std::size_t>(posL_) != opsL_.size() && !(opsL_[posL_] == dummyL)) ++posL_;
			
			//Test if pos = size 

			return std::log(std::abs(val_ = B_->at(posR_, posL_)));
		};
		
		int acceptErase() {
			int const N = opsL_.size(); int const newN = N - 1;
			
			if(newN) {
				int const inc = 1;
				Matrix* temp = new Matrix(newN, newN);
				
				if(posL_ != newN) { 
					dswap_(&N, B_->data(0, newN), &inc, B_->data(0, posL_), &inc); swapSign_ *= -1; 
					opsL_[posL_] = opsL_.back(); 					 
				}
				if(posR_ != newN) { 
					dswap_(&N, B_->data(newN, 0), &N, B_->data(posR_, 0), &N); swapSign_ *= -1;
					opsR_[posR_] = opsR_.back(); 
				} 
				
				for(int n = 0; n < newN; ++n) 
					dcopy_(&newN, B_->data(0, n), &inc, temp->data(0, n), &inc);
				
				double const fact = -1./val_;
				dger_(&newN, &newN, &fact, B_->data(0, newN), &inc, B_->data(newN, 0), &N, temp->data(), &newN);
				
				delete B_; B_ = temp;
			} else {
				delete B_; B_ = 0;
			}
			
			opsL_.pop_back(); opsR_.pop_back();
			
			det_ *= val_;
			
			return val_ > .0 ? 1 : -1;
		};
		
		template<class L>
		int rebuild(L const& link) {			
			int const N = opsL_.size();
			det_ = swapSign_;
			
			if(N) {
				if(!B_) B_ = new Matrix(N, N);
				
				Matrix toInvert(N, N);
				
				for(int j = 0; j < N; ++j) 					
					for(int i = 0; i < N; ++i) 
						toInvert.at(i,j) = link(opsL_[i], opsR_[j]);
				
				int const inc0 = 0; int const inc1 = 1;
				int const diagInc = N + 1; int const size = N*N;
				double const zero = .0; double const one = 1.;
				
				dcopy_(&size, &zero, &inc0, B_->data(), &inc1);
				dcopy_(&N, &one, &inc0, B_->data(), &diagInc);
				
				int ipiv[N]; int info;
				dgesv_(&N, &N, toInvert.data(), &N, ipiv, B_->data(), &N, &info);
				
				for(int i = 0; i < N; ++i) 
					det_ *= (ipiv[i] != i + 1 ? -toInvert.at(i, i) : toInvert.at(i, i));
			}
			
			return det_ > .0 ? 1 : -1;
		};
		
		double det() { return det_;};
		
		~Bath() { delete B_;};
	private:
		int swapSign_;
		Matrix* B_;
		double det_;
		
		std::vector<double> Bv_;
		std::vector<double> vec_;
		double val_;
		
		
		int posR_; 
		int posL_;
		Operator opR_; 
		Operator opL_;
		std::vector<Operator> opsR_;
		std::vector<Operator> opsL_;
	};
	
};

#endif






















