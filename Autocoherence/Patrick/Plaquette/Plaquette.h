#ifndef __PLAQUETTE
#define __PLAQUETTE

#include <boost/lexical_cast.hpp>
#include "../Flavors.h"

Fl::FlavorNames CuONames[] = {"d", "px", "py"};
typedef Fl::FlavorMatrix<3, CuONames, CuONames> CuOMatrix;

Fl::FlavorNames RCuNames[] = {"d_0", "d_1", "d_2", "d_3"};
typedef Fl::FlavorMatrix<4, RCuNames, RCuNames> RCuMatrix;

Fl::FlavorNames RCuONames[] = {"d_0", "px_0", "py_0", "d_1", "px_1", "py_1", "d_2", "px_2", "py_2", "d_3", "px_3", "py_3"};
typedef Fl::FlavorMatrix<12, RCuONames, RCuONames> RCuOMatrix;

Fl::FlavorNames RONames[] = {"px_0", "py_0", "px_1", "py_1", "px_2", "py_2", "px_3", "py_3"};
typedef Fl::FlavorMatrix<8, RONames, RONames> ROMatrix;

struct G0RCuInv {
	//Class for the non-interacting green's function, on the Cu site only
	G0RCuInv(double tpd, double tpp, double tppp, double ep) : tpd_(tpd), tpp_(tpp),tppp_(tppp), ep_(ep) {};
	
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
	double const tppp_;
	double const ep_;
	
	RCuMatrix result_;
	
	mutable CuOMatrix G0FullInv_;
	mutable CuOMatrix G0Full_;
	
	void add(std::complex<double> z, double kx, double ky, RCuMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		G0FullInv_        = CuOMatrix::Diag(z);  G0FullInv_(0, 1) -= tpd_*(1. - exp_mkx);            	G0FullInv_(0, 2) -= tpd_*(1. - exp_mky);
		G0FullInv_(1, 0) -= tpd_*(1. - exp_kx);  G0FullInv_(1, 1) -= ep_ - 2.*tpp_ + 2*tppp_*coskx;  	G0FullInv_(1, 2) -= tpp_*(1. - exp_kx)*(1. - exp_mky);
		G0FullInv_(2, 0) -= tpd_*(1. - exp_ky);  G0FullInv_(2, 1) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);	G0FullInv_(2, 2) -= ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
		
		G0FullInv_.inv(G0Full_);
		
		std::complex<double> v[] = {1., exp_kx, exp_kx*exp_ky, exp_ky};
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j)
				arg(i, j) += v[i]*std::conj(v[j])/G0Full_(0, 0);
	};
};

struct G0RCuOInv {
	//Class for the non-interacting green's function, on the Cu and O sites
	G0RCuOInv(double tpd, double tpp, double tppp, double ep) : tpd_(tpd), tpp_(tpp), tppp_(tppp), ep_(ep) {};
	
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
	double const tppp_;
	double const ep_;
	
	RCuOMatrix result_;
	
	mutable CuOMatrix G0FullInv_;
	
	void add(std::complex<double> z, double kx, double ky, RCuOMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		G0FullInv_        = CuOMatrix::Diag(z);  G0FullInv_(0, 1) -= tpd_*(1. - exp_mkx);                	G0FullInv_(0, 2) -= tpd_*(1. - exp_mky);
		G0FullInv_(1, 0) -= tpd_*(1. - exp_kx);  G0FullInv_(1, 1) -= ep_ - 2.*tpp_ + 2*tppp_*coskx;			G0FullInv_(1, 2) -= tpp_*(1. - exp_kx)*(1. - exp_mky);
		G0FullInv_(2, 0) -= tpd_*(1. - exp_ky);  G0FullInv_(2, 1) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);  	G0FullInv_(2, 2) -= ep_ - 2*tpp_ + 2.*tppp_*cosky;
		

		std::complex<double> v[] = {1., exp_kx, exp_kx*exp_ky, exp_ky};
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j) {
				std::complex<double> const fact = v[i]*std::conj(v[j]);
				for(int oi = 0; oi < 3; ++oi) 
					for(int oj = 0; oj < 3; ++oj)
				        arg(3*i + oi, 3*j + oj) += fact*G0FullInv_(oi, oj);
			};
	};
};

struct RCuLatticeGreen {
	//Class for the interacting Green's functzion (including the self energy in the computation)
	RCuLatticeGreen(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : z_(z), g0RCuInv_(tpd, tpp, tppp, ep), selfEnergy_(selfEnergy) {};
	
	RCuMatrix const& operator()(double kx, double ky) {
		temp_ = g0RCuInv_(z_, kx, ky);
		temp_ -= selfEnergy_;
		temp_.inv(result_);
		return result_;
	};
private:
	std::complex<double> const z_;
	G0RCuInv g0RCuInv_;
	RCuMatrix selfEnergy_;
	
	RCuMatrix temp_;
	RCuMatrix result_;
};

struct RCuOLatticeGreen {
	RCuOLatticeGreen(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : z_(z), g0RCuOInv_(tpd, tpp, tppp, ep), selfEnergy_(selfEnergy) {};
	
	RCuOMatrix const& operator()(double kx, double ky) {
		temp_ = g0RCuOInv_(z_, kx, ky);
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j)
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

struct RCuOLatticeHoppingMatrix {
	RCuOLatticeHoppingMatrix(double tpd, double tpp, double tppp, double ep) : tpd_(tpd), tpp_(tpp), tppp_(tppp), ep_(ep) {}; 
	
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
	double const tppp_;
	double const ep_;
	
	RCuOMatrix result_;
	
	mutable CuOMatrix t_;
	
	void add(double kx, double ky, RCuOMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
										t_(0, 1) = tpd_*(1. - exp_mkx);                t_(0, 2) = tpd_*(1. - exp_mky);
		t_(1, 0) = tpd_*(1. - exp_kx);  t_(1, 1) = ep_ - 2.*tpp_ + 2*tppp_*coskx;      t_(1, 2) = tpp_*(1. - exp_kx)*(1. - exp_mky); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		t_(2, 0) = tpd_*(1. - exp_ky);  t_(2, 1) = tpp_*(1. - exp_mkx)*(1. - exp_ky);  t_(2, 2) = ep_ - 2.*tpp_ + 2*tppp_*cosky;
		
		std::complex<double> v[] = {1., exp_kx, exp_kx*exp_ky, exp_ky};
		
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < 4; ++j) {
				std::complex<double> const fact = v[i]*std::conj(v[j]);
				for(int oi = 0; oi < 3; ++oi) 
					for(int oj = 0; oj < 3; ++oj)
				        arg(3*i + oi, 3*j + oj) += fact*t_(oi, oj);
			};
	};
};

struct RCuOLatticeKineticEnergy {
	RCuOLatticeKineticEnergy(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : latticeGreenRCuO_(z, tpd, tpp, tppp, ep, selfEnergy), latticeHoppingMatrixRCuO_(tpd, tpp, tppp, ep) {};
	
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