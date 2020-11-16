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
		Operator(int site, double time, int* ptr) : site_(site), time_(time), ptr_(ptr) {};
		int site() const { return site_;};
		double time() const { return time_;};
		int n() const { return *ptr_;};
	private:
		int site_;
		double time_;
		int* ptr_;
	};
	
	struct Bath {
		/** 
		* 
		* struct GreenIterator
		* 
		* Parameters :	opBegin : the operator which signifies the column of the Green function matrix
		*				opEnd : the end of the row 
		*				opDagg : the dagger operator (the rwo of the green's function matrix
		*				value : pointer on the value of the matrix at Green[opBegin,opDagg]
		* 
		* Description :
		*	This is a custom iterator over a Green's function matrix.
		*	This allow iteration over the Green function matrix in C++ standard while being able to know which component we are looking at (which operators are involved)
		*	opEnd allows the program to know when to go to the next row when reaching the last operators of the list to avoid Segmentation error or Access Violation.
		* 
		*/
		struct GreenIterator {
			GreenIterator(Operator const* opBegin, Operator const* opEnd, Operator const* opDagg, double const* value) : opBegin_(opBegin), opEnd_(opEnd), op_(opBegin), opDagg_(opDagg), value_(value) {};
			Operator const& op() const { return *op_;};
			Operator const& opDagg() const { return *opDagg_;};			
			double value() const { return *value_;};
			
			GreenIterator& operator++() { ++value_; if(++op_ == opEnd_) { op_ = opBegin_; ++opDagg_;}; return *this;}; 
			bool operator!=(GreenIterator const& other) const { return value_ != other.value_;};
		private:
			Operator const* const opBegin_; Operator const* const opEnd_;
			Operator const* op_; Operator const* opDagg_; double const* value_;
		};
		
		Bath() : swapSign_(1), B_(0), det_(1.) {};
		
		void add(int site, double time, int* ptr) { ops_.push_back(Operator(site, time, ptr));};
		void addDagg(int site, double time, int* ptr) { opsDagg_.push_back(Operator(site, time, ptr));};
		
		GreenIterator begin() const { return B_ ? GreenIterator(ops_.data(), ops_.data() + ops_.size(), opsDagg_.data(), B_->data()) : GreenIterator(0, 0, 0, 0);};
		GreenIterator end() const { return B_ ? GreenIterator(0, 0, 0, B_->data() + ops_.size()*opsDagg_.size()) : GreenIterator(0, 0, 0, 0);};
		
		/** 
		* 
		* template<typename L> double insert(int site, double time, int* ptr, double timeDagg, int* ptrDagg, L const& link)
		* 
		* Parameters :	site : site on which the vertex is added
		*				time : time of the annhilation operator
		*				ptr : pointer to the annihilation operator (only used to create the Operator)
		*				timeDagg : time of the creation operator
		*				ptrDagg :  pointer to the creation operator (only used to create the Operator)
		*				link : link function to get the hybridization function between sites
		* 
		* Return value : Returns log(|G[site,site]|) in the new configuration
		*
		* Description: 
		*   Computes G[site,site] in the new configuration if we added a site, dosen't affect the Green's matrix
		* 
		*/
		
		template<typename L> double insert(int site, double time, int* ptr, double timeDagg, int* ptrDagg, L const& link) {
			op_ = Operator(site, time, ptr); 
			opDagg_ = Operator(site, timeDagg, ptrDagg);
			
			int const N = opsDagg_.size();
			val_ = link(opDagg_.site(), op_.site(), opDagg_.time() - op_.time()); //Hybridization energy between the same site during the time the particle exists
			
			if(N) {  
				Bv_.resize(N); vec_.resize(N);
				for(int n = 0; n < N; ++n) vec_[n] = link(opsDagg_[n].site(), op_.site(), opsDagg_[n].time() - op_.time());
				
				char const no = 'n';
				int const inc = 1;
				double const zero = .0;
				double const one = 1.;
				dgemv_(&no, &N, &N, &one, B_->data(), &N, vec_.data(), &inc, &zero, Bv_.data(), &inc); // Bv_ = B_*vec_
				
				for(int n = 0; n < N; ++n) vec_[n] = link(opDagg_.site(), ops_[n].site(), opDagg_.time() - ops_[n].time());			
				val_ -= ddot_(&N, vec_.data(), &inc, Bv_.data(), &inc);	//Remove the hybridization energies of 		site->siteInMatrix->bath->siteInMatrix'->site 
			}
			
			return std::log(std::abs(val_));
		}; 
		


		/* 
		* int acceptInsert()
		* 
		* Return value : The sign of the G[site,site] coefficient
		*
		* Description: 
		*  	This function accepts the insertion of a new vertex. Therefore, we have to change the matrix B_ to include the new vertex. (the Green function Matrix)
		*	We enlarge the matrix, take into account the new vertex in the interactions between the old vertices (B_ + Bv_*(hBTilde)' on its N*N reduction)
		*	Then we add the interactions between the new vertex and the other vertices directly.
		*/
		int acceptInsert() {  
			int const N = opsDagg_.size();
			int const newN = N + 1;
			
			ops_.push_back(op_); opsDagg_.push_back(opDagg_); // on insère des  nouveaux operateurs de création et de destruction dans l'ensemble des operateurs
			
			Matrix* temp = new Matrix(newN, newN);
			
			double fact = 1./val_;
			temp->at(N, N) = fact;
			
			if(N) {                                    
				std::vector<double> hBTilde(N);
				
				char const yes = 't';
				int const inc = 1;
				double const zero = .0;
				double const one = 1.;
				dgemv_(&yes, &N, &N, &fact, B_->data(), &N, vec_.data(), &inc, &zero, hBTilde.data(), &inc);	//hBTilde = fact*(B_)'*vec_
				dger_(&N, &N, &one, Bv_.data(), &inc, hBTilde.data(), &inc, B_->data(), &N);					//B_ = B_ + Bv_*(hBTilde)'
				
				for(int n = 0; n < N; ++n) 
					dcopy_(&N, B_->data(0, n), &inc, temp->data(0, n), &inc);									//temp[:,n] = B_[:,n]
				
				fact = -1./val_;
				dscal_(&N, &fact, Bv_.data(), &inc); 															//Bv_*=fact
				dcopy_(&N, Bv_.data(), &inc, temp->data(0, N), &inc);											//temp[:,N] = Bv_
				
				double const minus = -1.;
				dscal_(&N, &minus, hBTilde.data(), &inc); 														//hBTilde*=-1
				dcopy_(&N, hBTilde.data(), &inc, temp->data(N, 0), &newN); 										//temp[N,:] = hBTilde
			}
			
			delete B_; B_ = temp;
			
			det_ *= val_;
			
			return val_ > .0 ? 1 : -1;
		};
		/** 
		* 
		* double erase(int site, double time, double timeDagg)
		* 
		* Parameters :	site : Site of the vertex to remove
						time : time of the annhilation operator
						timeDagg : time of the creation operator
		* 
		* Return value : Returns log(|G[site,site]|) in the old configuration
		*
		* Description: 
		*   Prepares the erase by getting the right vertex(using site and time to look for it) and the value of it's green function
		*	We keep in memory which vertex we were tryuing to erase
		* 
		*/
		double erase(int site, double time, double timeDagg) {
			pos_ = 0; while(static_cast<std::size_t>(pos_) != ops_.size() && !(ops_[pos_].time() == time && ops_[pos_].site() == site)) ++pos_;
			posDagg_ = 0; while(static_cast<std::size_t>(posDagg_) != opsDagg_.size() && !(opsDagg_[posDagg_].time() == timeDagg && opsDagg_[posDagg_].site() == site)) ++posDagg_;

			return std::log(std::abs(val_ = B_->at(pos_, posDagg_)));
		};
		
		/* 
		* int acceptErase()
		* 
		* Return value : The sign of the old G[site,site] coefficient
		*
		* Description: 
		*  	This function accepts the erase of a new vertex. Therefore, we have to change the matrix B_ to remove this vertex. (the Green function Matrix)
		*	Here, we want to make B_ smaller.
		*	We then substract the action of the old vertex on the system in this smaller matrix and update the determinant of the matrix.
		*/
		int acceptErase() {
			int const N = opsDagg_.size(); int const newN = N - 1;
			
			if(newN) {
				int const inc = 1;
				Matrix* temp = new Matrix(newN, newN);
				
				//If need be, we exchange the rows in the B_ matrix to put the informations relative to posDagg_ and pos_ on the newN+1 rows and columns, to get rid of them more easily
				if(posDagg_ != newN) { 
					dswap_(&N, B_->data(0, newN), &inc, B_->data(0, posDagg_), &inc); swapSign_ *= -1; 				
					opsDagg_[posDagg_] = opsDagg_.back(); 					 
				}
				if(pos_ != newN) { 
					dswap_(&N, B_->data(newN, 0), &N, B_->data(pos_, 0), &N); swapSign_ *= -1;
					ops_[pos_] = ops_.back(); 
				} 
				
				//We get rid of the last row and the last column of B_
				for(int n = 0; n < newN; ++n) 
					dcopy_(&newN, B_->data(0, n), &inc, temp->data(0, n), &inc);
				
				//For some reason, we also have to account for the removed position in the new B
				double const fact = -1./val_;
				dger_(&newN, &newN, &fact, B_->data(0, newN), &inc, B_->data(newN, 0), &N, temp->data(), &newN);    // temp = temp + fact*B_[*,newN]*B_[newN,*]
				
				delete B_; B_ = temp;
			} else {
				delete B_; B_ = 0;
			}
			
			opsDagg_.pop_back(); ops_.pop_back();
			
			det_ *= val_;
			
			return val_ > .0 ? 1 : -1;
		};
		
		/* 
		* template<class L>
		* int rebuild(L const& link)
		* 
		* Return value : The sign of the matrix determinant
		*
		* Description: 
		*  	This function rebuilds the Green function Matrix from the link function. It also computes the determinant of the matrix
		*/
		template<class L>
		int rebuild(L const& link) {			
			int const N = opsDagg_.size();
			det_ = swapSign_;
			
			if(N) {
				if(!B_) B_ = new Matrix(N, N);
				
				Matrix toInvert(N, N);
				
				for(int j = 0; j < N; ++j) 					
					for(int i = 0; i < N; ++i) 
						toInvert.at(i,j) = link(opsDagg_[i].site(), ops_[j].site(), opsDagg_[i].time() - ops_[j].time());
				
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
		
		
		int pos_; 
		int posDagg_;
		Operator op_; 
		Operator opDagg_;
		std::vector<Operator> ops_;
		std::vector<Operator> opsDagg_;
	};
	
};

#endif






















