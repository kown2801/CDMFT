#ifndef __PLAQUETTE
#define __PLAQUETTE

#include <boost/lexical_cast.hpp>
#include "../Flavors.h"

Fl::FlavorNames RCuNames[] = {"d_0Up", "d_1Up", "d_2Up", "d_3Up", "d_0Down", "d_1Down", "d_2Down", "d_3Down"};
Fl::FlavorNames CuONames[] = {"d", "px", "py"};

typedef Fl::FlavorMatrix<8, RCuNames, RCuNames> RCuMatrix;
typedef Fl::FlavorMatrix<3, CuONames, CuONames> CuOMatrix;

struct G0RCuInv {
	G0RCuInv(double tpd, double tpp, double ep) : tpd_(tpd), tpp_(tpp), ep_(ep) {};
	
	RCuMatrix const& operator()(std::complex<double> z, double kx, double ky) {
		result_ = .0;
		add(z, kx, ky, result_);
		add(z, kx + M_PI, ky, result_);
		add(z, kx + M_PI, ky + M_PI, result_);
		add(z, kx, ky + M_PI, result_);
		result_ *= .25;
		return result_;
	};
private:
	double const tpd_;
	double const tpp_;
	double const ep_;
	
	RCuMatrix result_;
	
	mutable CuOMatrix G0kFullInv_;
	mutable CuOMatrix G0mkFullInv_;
	mutable CuOMatrix G0kFull_;
	mutable CuOMatrix G0mkFull_;
	
	void add(std::complex<double> z, double kx, double ky, RCuMatrix& arg) {		
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		G0kFullInv_        = CuOMatrix::Diag(z);  G0kFullInv_(0, 1) -= tpd_*(1. - exp_mkx);                G0kFullInv_(0, 2) -= tpd_*(1. - exp_mky);
		G0kFullInv_(1, 0) -= tpd_*(1. - exp_kx);  G0kFullInv_(1, 1) -= ep_ + 2.*tpp_*(coskx - 1.);         G0kFullInv_(1, 2) -= tpp_*(1. - exp_kx)*(1. - exp_mky);
		G0kFullInv_(2, 0) -= tpd_*(1. - exp_ky);  G0kFullInv_(2, 1) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);  G0kFullInv_(2, 2) -= ep_ + 2.*tpp_*(cosky - 1.);
		
		G0mkFullInv_        = CuOMatrix::Diag(z);   G0mkFullInv_(0, 1) -= tpd_*(1. - exp_kx);                 G0mkFullInv_(0, 2) -= tpd_*(1. - exp_ky);
		G0mkFullInv_(1, 0) -= tpd_*(1. - exp_mkx);  G0mkFullInv_(1, 1) -= ep_ + 2.*tpp_*(coskx - 1.);         G0mkFullInv_(1, 2) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);
		G0mkFullInv_(2, 0) -= tpd_*(1. - exp_mky);  G0mkFullInv_(2, 1) -= tpp_*(1. - exp_kx)*(1. - exp_mky);  G0mkFullInv_(2, 2) -= ep_ + 2.*tpp_*(cosky - 1.);
		
		G0kFullInv_.inv(G0kFull_);
		
		G0mkFullInv_.inv(G0mkFull_);
		
		std::complex<double> v[] = {1., exp_kx, exp_kx*exp_ky, exp_ky};
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j) {
				arg(i, j) += v[i]*std::conj(v[j])/G0kFull_(0, 0);
				arg(i + 4, j + 4) += v[i]*std::conj(v[j])/(-std::conj(G0mkFull_(0, 0)));  
			};
	};
};


Fl::FlavorNames RCuONames[] = {"d_0Up", "px_0Up", "py_0Up", "d_1Up", "px_1Up", "py_1Up", "d_2Up", "px_2Up", "py_2Up", "d_3Up", "px_3Up", "py_3Up",
                               "d_0Down", "px_0Down", "py_0Down", "d_1Down", "px_1Down", "py_1Down", "d_2Down", "px_2Down", "py_2Down", "d_3Down", "px_3Down", "py_3Down"};

typedef Fl::FlavorMatrix<24, RCuONames, RCuONames> RCuOMatrix;


struct G0RCuOInv {
	G0RCuOInv(double tpd, double tpp, double ep) : tpd_(tpd), tpp_(tpp), ep_(ep) {};
	
	RCuOMatrix const& operator()(std::complex<double> z, double kx, double ky) {
		result_ = .0;
		add(z, kx, ky, result_);
		add(z, kx + M_PI, ky, result_);
		add(z, kx + M_PI, ky + M_PI, result_);
		add(z, kx, ky + M_PI, result_);
		result_ *= .25;
		return result_;
	};
private:
	double const tpd_;
	double const tpp_;
	double const ep_;
	
	RCuOMatrix result_;

	mutable CuOMatrix G0kFullInv_;
	mutable CuOMatrix G0mkFullInv_;
	
	void add(std::complex<double> z, double kx, double ky, RCuOMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		G0kFullInv_        = CuOMatrix::Diag(z);  G0kFullInv_(0, 1) -= tpd_*(1. - exp_mkx);                G0kFullInv_(0, 2) -= tpd_*(1. - exp_mky);
		G0kFullInv_(1, 0) -= tpd_*(1. - exp_kx);  G0kFullInv_(1, 1) -= ep_ + 2.*tpp_*(coskx - 1.);         G0kFullInv_(1, 2) -= tpp_*(1. - exp_kx)*(1. - exp_mky);
		G0kFullInv_(2, 0) -= tpd_*(1. - exp_ky);  G0kFullInv_(2, 1) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);  G0kFullInv_(2, 2) -= ep_ + 2.*tpp_*(cosky - 1.);
		
		G0mkFullInv_        = CuOMatrix::Diag(z);   G0mkFullInv_(0, 1) -= tpd_*(1. - exp_kx);                 G0mkFullInv_(0, 2) -= tpd_*(1. - exp_ky);
		G0mkFullInv_(1, 0) -= tpd_*(1. - exp_mkx);  G0mkFullInv_(1, 1) -= ep_ + 2.*tpp_*(coskx - 1.);         G0mkFullInv_(1, 2) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);
		G0mkFullInv_(2, 0) -= tpd_*(1. - exp_mky);  G0mkFullInv_(2, 1) -= tpp_*(1. - exp_kx)*(1. - exp_mky);  G0mkFullInv_(2, 2) -= ep_ + 2.*tpp_*(cosky - 1.);
		
		std::complex<double> v[] = {1., exp_kx, exp_kx*exp_ky, exp_ky};
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j) {
				std::complex<double> const fact = v[i]*std::conj(v[j]);
				for(int oi = 0; oi < 3; ++oi) 
					for(int oj = 0; oj < 3; ++oj) {
				        arg(3*i + oi, 3*j + oj) += fact*G0kFullInv_(oi, oj);
						
						arg(3*(i + 4) + oi, 3*(j + 4) + oj) += fact*(-std::conj(G0mkFullInv_(oi, oj)));
					}
			};		
	};
};


struct RCuLatticeGreen {
	RCuLatticeGreen(std::complex<double> z, double tpd, double tpp, double ep, RCuMatrix const& selfEnergy) : z_(z), g0RInv_(tpd, tpp, ep), selfEnergy_(selfEnergy) {};
	RCuMatrix const& operator()(double kx, double ky) {
		temp_ = g0RInv_(z_, kx, ky);
		temp_ -= selfEnergy_;
		temp_.inv(result_);
		return result_;
	};
private:
	std::complex<double> const z_;
	G0RCuInv g0RInv_;
	RCuMatrix selfEnergy_;
	
	RCuMatrix temp_;
	RCuMatrix result_;
};

struct RCuOLatticeGreen {
	RCuOLatticeGreen(std::complex<double> z, double tpd, double tpp, double ep, RCuMatrix const& selfEnergy) : z_(z), g0RCuOInv_(tpd, tpp, ep), selfEnergy_(selfEnergy) {};
	RCuOMatrix const& operator()(double kx, double ky) {
		temp_ = g0RCuOInv_(z_, kx, ky);
		
		for(int i = 0; i < 8; ++i)
			for(int j = 0; j < 8; ++j)
				temp_(3*i, 3*j) -= selfEnergy_(i, j);
		
		temp_.inv(result_);
		return result_;
	};
private:
	std::complex<double> const z_;
	G0RCuOInv g0RCuOInv_;
	RCuMatrix selfEnergy_;
	
	RCuOMatrix temp_;
	RCuOMatrix result_;
};

Fl::FlavorNames RONames[] = {"px_0Up", "py_0Up", "px_1Up", "py_1Up", "px_2Up", "py_2Up", "px_3Up", "py_3Up", 
                             "px_0Down", "py_0Down", "px_1Down", "py_1Down", "px_2Down", "py_2Down", "px_3Down", "py_3Down"};

typedef Fl::FlavorMatrix<16, RONames, RONames> ROMatrix;


struct RCuOLatticeHoppingMatrix {
	RCuOLatticeHoppingMatrix(double tpd, double tpp, double ep) : tpd_(tpd), tpp_(tpp), ep_(ep) {}; 
	
	RCuOMatrix const& operator()(double kx, double ky) {
		result_ = .0;
		add(kx, ky, result_);
		add(kx + M_PI, ky, result_);
		add(kx + M_PI, ky + M_PI, result_);
		add(kx, ky + M_PI, result_);
		result_ *= .25;
		
		return result_;
	};
	
private:
	double const tpd_;
	double const tpp_;
	double const ep_;
	
	RCuOMatrix result_;
	
	mutable CuOMatrix tk_;
	mutable CuOMatrix tmk_;
	
	void add(double kx, double ky, RCuOMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		                                 tk_(0, 1) = tpd_*(1. - exp_mkx);                tk_(0, 2) = tpd_*(1. - exp_mky);
		tk_(1, 0) = tpd_*(1. - exp_kx);  tk_(1, 1) = ep_ + 2.*tpp_*(coskx - 1.);         tk_(1, 2) = tpp_*(1. - exp_kx)*(1. - exp_mky);
		tk_(2, 0) = tpd_*(1. - exp_ky);  tk_(2, 1) = tpp_*(1. - exp_mkx)*(1. - exp_ky);  tk_(2, 2) = ep_ + 2.*tpp_*(cosky - 1.);
		
		                                   tmk_(0, 1) = tpd_*(1. - exp_kx);                 tmk_(0, 2) = tpd_*(1. - exp_ky);
		tmk_(1, 0) = tpd_*(1. - exp_mkx);  tmk_(1, 1) = ep_ + 2.*tpp_*(coskx - 1.);         tmk_(1, 2) = tpp_*(1. - exp_mkx)*(1. - exp_ky);
		tmk_(2, 0) = tpd_*(1. - exp_mky);  tmk_(2, 1) = tpp_*(1. - exp_kx)*(1. - exp_mky);  tmk_(2, 2) = ep_ + 2.*tpp_*(cosky - 1.);
		
		std::complex<double> v[] = {1., exp_kx, exp_kx*exp_ky, exp_ky};
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j) {
				std::complex<double> const fact = v[i]*std::conj(v[j]);
				for(int oi = 0; oi < 3; ++oi) 
					for(int oj = 0; oj < 3; ++oj) {
				        arg(3*i + oi, 3*j + oj) += fact*tk_(oi, oj);
						
						arg(3*(i + 4) + oi, 3*(j + 4) + oj) += fact*(-std::conj(tmk_(oi, oj)));
					}
			};
	};
};

struct RCuOLatticeKineticEnergy {
	RCuOLatticeKineticEnergy(std::complex<double> z, double tpd, double tpp, double ep, RCuMatrix const& selfEnergy) : latticeGreenRCuO_(z, tpd, tpp, ep, selfEnergy), latticeHoppingMatrixRCuO_(tpd, tpp, ep) {};
	
	RCuOMatrix const& operator()(double kx, double ky) {
        result_ = latticeHoppingMatrixRCuO_(kx, ky)*latticeGreenRCuO_(kx, ky);
		
		return result_;
	};
	
private:
	RCuOLatticeGreen latticeGreenRCuO_;
	RCuOLatticeHoppingMatrix latticeHoppingMatrixRCuO_;
	
	RCuOMatrix result_;
};

#endif





