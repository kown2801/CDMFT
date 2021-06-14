#ifndef __FLAVORS
#define __FLAVORS

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <set>
#include <bitset>
#include <algorithm>
#include <fstream>
#include <string>

extern "C" void   zgesv_(const int*, const int*, std::complex<double>*, const int*, int*, std::complex<double>*, const int*, int*);

namespace Fl {
	typedef const char* FlavorNames;
	
	typedef unsigned long State;
	
	template<int N, FlavorNames (&names)[N]>
	struct FlavorState : public std::bitset<N> {
		FlavorState() : sign_(1) {};
		explicit FlavorState(State const& state) : std::bitset<N>(state), sign_(1) {};
		FlavorState(FlavorState const& state) : std::bitset<N>(state), sign_(state.sign()) {};
		
		State state() const { return this->to_ulong();};
		int sign() const { return sign_;};
		int& sign() { return sign_;};
		
		int operator()(std::string str) const { 
			const char** it = std::find(names, names + N, str);
			if(it == names + N) throw std::runtime_error("FlavorState: flavor " + str + " not found.");
			return this->test(it - names);
		};
	private:
		int sign_;
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N]> struct FlavorMatrix;
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N]>
	struct Hc {
		Hc(FlavorMatrix<N, namesI, namesJ> const& source) : source_(source) {};
		std::complex<double> operator()(int i, int j) const { return std::conj(source_(j, i));};
		void equ(FlavorMatrix<N, namesJ, namesI>& dest) const { for(int i = 0; i < N; ++i) for(int j = 0; j < N; ++j) dest(i, j) = std::conj(source_(j, i));};
		void mequ(FlavorMatrix<N, namesJ, namesI>& dest) const { for(int i = 0; i < N; ++i) for(int j = 0; j < N; ++j) dest(i, j) -= std::conj(source_(j, i));};
		void pequ(FlavorMatrix<N, namesJ, namesI>& dest) const { for(int i = 0; i < N; ++i) for(int j = 0; j < N; ++j) dest(i, j) += std::conj(source_(j, i));};
	private:
		FlavorMatrix<N, namesI, namesJ> const& source_;
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N], class Op>
	struct Scal {
		Scal(std::complex<double> val, Op const& source) : val_(val), source_(source) {};
		void equ(FlavorMatrix<N, namesI, namesJ>& dest) const { for(int j = 0; j < N; ++j) for(int i = 0; i < N; ++i) dest(i, j) = val_*source_(i, j);};
		void mequ(FlavorMatrix<N, namesI, namesJ>& dest) const { for(int j = 0; j < N; ++j) for(int i = 0; i < N; ++i) dest(i, j) -= val_*source_(i, j);};
		void pequ(FlavorMatrix<N, namesI, namesJ>& dest) const { for(int j = 0; j < N; ++j) for(int i = 0; i < N; ++i) dest(i, j) += val_*source_(i, j);};
	private:
		std::complex<double> const val_; 
		Op const& source_;
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N], class OpL, class OpR>
	struct Mult {
		Mult(OpL const& L, OpR const& R) : L_(L), R_(R) {};
		void equ(FlavorMatrix<N, namesI, namesJ>& dest) const { for(int i = 0; i < N; ++i) for(int j = 0; j < N; ++j) { std::complex<double> temp(.0); for(int k = 0; k < N; ++k) temp += L_(i, k)*R_(k, j); dest(i, j) = temp;};};
		void mequ(FlavorMatrix<N, namesI, namesJ>& dest) const { for(int i = 0; i < N; ++i) for(int j = 0; j < N; ++j) { std::complex<double> temp(.0); for(int k = 0; k < N; ++k) temp += L_(i, k)*R_(k, j); dest(i, j) -= temp;};};
		void pequ(FlavorMatrix<N, namesI, namesJ>& dest) const { for(int i = 0; i < N; ++i) for(int j = 0; j < N; ++j) { std::complex<double> temp(.0); for(int k = 0; k < N; ++k) temp += L_(i, k)*R_(k, j); dest(i, j) += temp;};};
	private:
		OpL const& L_;
		OpR const& R_;
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N]>
	struct FlavorMatrix {
		struct Diag {
			Diag(std::complex<double> val) : val_(val) {};
			void equ(FlavorMatrix& dest) const { for(int i = 0; i < N; ++i) dest(i, i) = val_;};
			void mequ(FlavorMatrix& dest) const { for(int i = 0; i < N; ++i) dest(i, i) -= val_;};
			void pequ(FlavorMatrix& dest) const { for(int i = 0; i < N; ++i) dest(i, i) += val_;};
		private:
			std::complex<double> const val_;
		};
		
		//-------------------------------------------------------------------------------------------------------------------------
		FlavorMatrix() { *this = .0;};
		FlavorMatrix(int val) { *this = val;};
		FlavorMatrix(double val) { *this = val;};
		FlavorMatrix(std::complex<double> val) { *this = val;};
		FlavorMatrix(FlavorMatrix const& matrix) { matrix.equ(*this);};
		
		FlavorMatrix(Diag const& diag) { *this = .0; diag.equ(*this);};
		template<class T> FlavorMatrix(T const& t) { t.equ(*this);};
		
		std::complex<double>& operator()(int i, int j) { return data_[i + N*j];};
		std::complex<double> const& operator()(int i, int j) const { return data_[i + N*j];};
		std::complex<double>& operator()(std::string strI, std::string strJ) { return data_[getIndexI(strI) + N*getIndexJ(strJ)];};
		std::complex<double> const& operator()(std::string strI, std::string strJ) const { return data_[getIndexI(strI) + N*getIndexJ(strJ)];};
		
		FlavorMatrix& operator=(int val) { for(int i = 0; i < N*N; ++i) data_[i] = val; return *this;};
		FlavorMatrix& operator=(double val) { for(int i = 0; i < N*N; ++i) data_[i] = val; return *this;};
		FlavorMatrix& operator=(std::complex<double> val) { for(int i = 0; i < N*N; ++i) data_[i] = val; return *this;};
		FlavorMatrix& operator*=(std::complex<double> fact) { for(int i = 0; i < N*N; ++i) data_[i] *= fact; return *this;};
		FlavorMatrix& operator/=(std::complex<double> divisor) { for(int i = 0; i < N*N; ++i) data_[i] /= divisor; return *this;};
		
		void equ(FlavorMatrix& dest) const { for(int i = 0; i < N*N; ++i) dest.data_[i] = data_[i];}; 
		void mequ(FlavorMatrix& dest) const { for(int i = 0; i < N*N; ++i) dest.data_[i] -= data_[i];}; 
		void pequ(FlavorMatrix& dest) const { for(int i = 0; i < N*N; ++i) dest.data_[i] += data_[i];};
		
		Hc<N, namesI, namesJ> hc() const { return Hc<N, namesI, namesJ>(*this);};
		
		int size(){return N;};
		
		FlavorMatrix& operator=(Diag const& diag) { *this = .0; diag.equ(*this); return *this;};		
		template<class T> FlavorMatrix& operator=(T const& t) { t.equ(*this); return *this;};
		template<class T> FlavorMatrix& operator-=(T const& t) { t.mequ(*this); return *this;};
		template<class T> FlavorMatrix& operator+=(T const& t) { t.pequ(*this); return *this;};
		
		FlavorMatrix& inv() {
			FlavorMatrix temp = *this; *this = Diag(1.);
			
			int ipiv[N]; int info; int dim = N;
			zgesv_(&dim, &dim, temp.data_, &dim, ipiv, data_, &dim, &info);
			
			return *this;
		};
		
		void inv(FlavorMatrix& result) {
			result = Diag(1.);
			
			int ipiv[N]; int info; int dim = N;
			zgesv_(&dim, &dim, data_, &dim, ipiv, result.data_, &dim, &info);
		};
		
		double abs() const {
			double temp = .0; for(int i = 0; i < N*N; ++i) temp += data_[i].real()*data_[i].real() + data_[i].imag()*data_[i].imag();
			return std::sqrt(temp);
		};
		
		std::complex<double> trace() const {
			std::complex<double> temp = .0; for(int i = 0; i < N; ++i) temp += data_[(N + 1)*i]; return temp;
		};
		
		~FlavorMatrix() {};
	private:
		std::complex<double> data_[N*N];
		
		int getIndexI(std::string const& str) const {
			const char** it = std::find(namesI, namesI + N, str);
			if(it == namesI + N) throw std::runtime_error("FlavorMatrix: flavor " + str + " not found.");
			return it - namesI;
		};
		
		int getIndexJ(std::string const& str) const {
			const char** it = std::find(namesJ, namesJ + N, str);
			if(it == namesJ + N) throw std::runtime_error("FlavorMatrix: flavor " + str + " not found.");
			return it - namesJ;
		};
		
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N]>
	Scal<N, namesI, namesJ, FlavorMatrix<N, namesI, namesJ> > operator*(std::complex<double> const& fact, FlavorMatrix<N, namesI, namesJ> const& matrix) {
		return Scal<N, namesI, namesJ, FlavorMatrix<N, namesI, namesJ> >(fact, matrix);
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N]>
	Scal<N, namesI, namesJ, Hc<N, namesJ, namesI> > operator*(std::complex<double> const& fact, Hc<N, namesJ, namesI> const& matrix) {
		return Scal<N, namesI, namesJ, Hc<N, namesJ, namesI> >(fact, matrix);
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesK)[N], FlavorNames (&namesJ)[N]>
	Mult<N, namesI, namesJ, FlavorMatrix<N, namesI, namesK>, FlavorMatrix<N, namesK, namesJ> > operator*(FlavorMatrix<N, namesI, namesK> const& L, FlavorMatrix<N, namesK, namesJ> const& R) {
		return Mult<N, namesI, namesJ, FlavorMatrix<N, namesI, namesK>, FlavorMatrix<N, namesK, namesJ> >(L, R);
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesK)[N], FlavorNames (&namesJ)[N]>
	Mult<N, namesI, namesJ, Hc<N, namesK, namesI>, FlavorMatrix<N, namesK, namesJ> > operator*(Hc<N, namesK, namesI> const& L, FlavorMatrix<N, namesK, namesJ> const& R) {
		return Mult<N, namesI, namesJ, Hc<N, namesK, namesI>, FlavorMatrix<N, namesK, namesJ> >(L, R);
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesK)[N], FlavorNames (&namesJ)[N]>
	Mult<N, namesI, namesJ, FlavorMatrix<N, namesI, namesK>, Hc<N, namesJ, namesK> > operator*(FlavorMatrix<N, namesI, namesK> const& L, Hc<N, namesJ, namesK> const& R) {
		return Mult<N, namesI, namesJ, FlavorMatrix<N, namesI, namesK>, Hc<N, namesJ, namesK> >(L, R);
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesK)[N], FlavorNames (&namesJ)[N]>
	Mult<N, namesI, namesJ, Hc<N, namesK, namesI>, Hc<N, namesJ, namesK> > operator*(Hc<N, namesK, namesI> const& L, Hc<N, namesJ, namesK> const& R) {
		return Mult<N, namesI, namesJ, Hc<N, namesK, namesI>, Hc<N, namesJ, namesK> >(L, R);
	};
	
	template<int N, FlavorNames (&namesI)[N], FlavorNames (&namesJ)[N]>
	double abs(FlavorMatrix<N, namesI, namesJ> const& arg) {
		return arg.abs();
	};	
};

#endif