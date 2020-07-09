#ifndef __PLAQUETTE
#define __PLAQUETTE

#include <boost/lexical_cast.hpp>
#include "../Flavors.h"

Fl::FlavorNames RCuNames[] = {"d_0Up", "d_1Up", "d_2Up", "d_3Up", "d_0Down", "d_1Down", "d_2Down", "d_3Down"};
Fl::FlavorNames CuONames[] = {"d", "px", "py"};
Fl::FlavorNames CuOSNames[] = {"d_Up","px_Up","py_Up","d_Down","px_Down","py_Down"};

typedef Fl::FlavorMatrix<8, RCuNames, RCuNames> RCuMatrix;
typedef Fl::FlavorMatrix<3, CuONames, CuONames> CuOMatrix;
typedef Fl::FlavorMatrix<6, CuOSNames, CuOSNames> CuOSMatrix;

struct G0RCuInv {
	G0RCuInv(double tpd, double tpp, double tppp, double ep) : tpd_(tpd), tpp_(tpp), tppp_(tppp), ep_(ep) {};
	
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
		
		G0kFullInv_        = CuOMatrix::Diag(z);  G0kFullInv_(0, 1) -= tpd_*(1. - exp_mkx);              G0kFullInv_(0, 2) -= tpd_*(1. - exp_mky);
		G0kFullInv_(1, 0) -= tpd_*(1. - exp_kx);  G0kFullInv_(1, 1) -= ep_ - 2*tpp_ + 2.*tppp_*coskx;     G0kFullInv_(1, 2) -= tpp_*(1. - exp_kx)*(1. - exp_mky);
		G0kFullInv_(2, 0) -= tpd_*(1. - exp_ky);  G0kFullInv_(2, 1) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);G0kFullInv_(2, 2) -= ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
		G0mkFullInv_        = CuOMatrix::Diag(z);   G0mkFullInv_(0, 1) -= tpd_*(1. - exp_kx);          		G0mkFullInv_(0, 2) -= tpd_*(1. - exp_ky);
		G0mkFullInv_(1, 0) -= tpd_*(1. - exp_mkx);  G0mkFullInv_(1, 1) -= ep_ - 2*tpp_ + 2.*tppp_*coskx;    	G0mkFullInv_(1, 2) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);
		G0mkFullInv_(2, 0) -= tpd_*(1. - exp_mky);  G0mkFullInv_(2, 1) -= tpp_*(1. - exp_kx)*(1. - exp_mky);G0mkFullInv_(2, 2) -= ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
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

	mutable CuOMatrix G0kFullInv_;
	mutable CuOMatrix G0mkFullInv_;
	
	void add(std::complex<double> z, double kx, double ky, RCuOMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		G0kFullInv_        = CuOMatrix::Diag(z);  G0kFullInv_(0, 1) -= tpd_*(1. - exp_mkx);             	G0kFullInv_(0, 2) -= tpd_*(1. - exp_mky);
		G0kFullInv_(1, 0) -= tpd_*(1. - exp_kx);  G0kFullInv_(1, 1) -= ep_ - 2*tpp_ + 2.*tppp_*coskx;    	G0kFullInv_(1, 2) -= tpp_*(1. - exp_kx)*(1. - exp_mky);
		G0kFullInv_(2, 0) -= tpd_*(1. - exp_ky);  G0kFullInv_(2, 1) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);  	G0kFullInv_(2, 2) -= ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
		G0mkFullInv_        = CuOMatrix::Diag(z);   G0mkFullInv_(0, 1) -= tpd_*(1. - exp_kx);               G0mkFullInv_(0, 2) -= tpd_*(1. - exp_ky);
		G0mkFullInv_(1, 0) -= tpd_*(1. - exp_mkx);  G0mkFullInv_(1, 1) -= ep_ - 2*tpp_ + 2.*tppp_*coskx;    	G0mkFullInv_(1, 2) -= tpp_*(1. - exp_mkx)*(1. - exp_ky);
		G0mkFullInv_(2, 0) -= tpd_*(1. - exp_mky);  G0mkFullInv_(2, 1) -= tpp_*(1. - exp_kx)*(1. - exp_mky);G0mkFullInv_(2, 2) -= ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
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
	RCuLatticeGreen(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : z_(z), g0RInv_(tpd, tpp, tppp, ep), selfEnergy_(selfEnergy) {};
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
	RCuOLatticeGreen(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : z_(z), g0RCuOInv_(tpd, tpp, tppp, ep), selfEnergy_(selfEnergy) {};
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
	
	mutable CuOMatrix tk_;
	mutable CuOMatrix tmk_;
	
	void add(double kx, double ky, RCuOMatrix& arg) {
		double coskx = std::cos(kx);
		std::complex<double> exp_kx(coskx, std::sin(kx));
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		std::complex<double> exp_ky(cosky, std::sin(ky));
		std::complex<double> exp_mky = std::conj(exp_ky);
		
		                                 tk_(0, 1) = tpd_*(1. - exp_mkx);               tk_(0, 2) = tpd_*(1. - exp_mky);
		tk_(1, 0) = tpd_*(1. - exp_kx);  tk_(1, 1) = ep_ - 2*tpp_ + 2.*tppp_*coskx;    	tk_(1, 2) = tpp_*(1. - exp_kx)*(1. - exp_mky);
		tk_(2, 0) = tpd_*(1. - exp_ky);  tk_(2, 1) = tpp_*(1. - exp_mkx)*(1. - exp_ky); tk_(2, 2) = ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
		                                   tmk_(0, 1) = tpd_*(1. - exp_kx);                 tmk_(0, 2) = tpd_*(1. - exp_ky);
		tmk_(1, 0) = tpd_*(1. - exp_mkx);  tmk_(1, 1) = ep_ - 2*tpp_ + 2.*tppp_*coskx;    	tmk_(1, 2) = tpp_*(1. - exp_mkx)*(1. - exp_ky);
		tmk_(2, 0) = tpd_*(1. - exp_mky);  tmk_(2, 1) = tpp_*(1. - exp_kx)*(1. - exp_mky);  tmk_(2, 2) = ep_ - 2*tpp_ + 2.*tppp_*cosky;
		
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

struct RCuOLatticeGreenPeriodized{
	RCuOLatticeGreenPeriodized(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : gSuperLattice_(z,tpd, tpp, tppp, ep,selfEnergy) {};
		CuOSMatrix& operator()(double kx, double ky) {
			std::complex<double> exp_kx(std::cos(kx), std::sin(kx));
			std::complex<double> exp_ky(std::cos(ky), std::sin(ky));
			if(kx > M_PI/2 || kx < -M_PI/2){
				kx-=M_PI;
			}
			if(ky > M_PI/2 || ky < -M_PI/2){
				ky-=M_PI;
			}
			temp_ = gSuperLattice_(kx,ky);
			result_ = 0;

			std::complex<double> v[] = {1., std::conj(exp_kx), std::conj(exp_kx*exp_ky), std::conj(exp_ky)};
			for(int i = 0; i < 4; ++i){
				for(int j = 0; j < 4; ++j){
					std::complex<double> const fact = v[i]*std::conj(v[j]);
					for(int oi = 0; oi < 3; ++oi){
						for(int oj = 0; oj < 3; ++oj) {
					        result_(oi,oj) += fact*temp_(3*i + oi, 3*j + oj);
					        result_(oi + 3,oj) += fact*temp_(3*(i+4) + oi, 3*j + oj);
					        result_(oi,oj + 3) += fact*temp_(3*i + oi, 3*(j+4) + oj);
					        result_(oi+ 3,oj + 3) += fact*temp_(3*(i+4) + oi, 3*(j+4) + oj);
						}
					}
				}
			}
			result_*=0.25;
			return result_;
		};
private:
	RCuOLatticeGreen gSuperLattice_;
	CuOSMatrix result_;

	RCuOMatrix temp_;
};

struct SuperfluidStiffness {
	SuperfluidStiffness(std::complex<double> z, double tpd, double tpp, double tppp, double ep, RCuMatrix const& selfEnergy) : tpd_(tpd), tpp_(tpp), tppp_(tppp), ep_(ep),z_(z),totalGreen_(z,tpd, tpp, tppp, ep,selfEnergy) {};
	
	std::complex<double> operator()(double kx, double ky) {
		CuOSMatrix green = totalGreen_(kx,ky);
		CuOSMatrix first_d_x = construct_first_d_x(kx,ky);
		CuOSMatrix tau3 = construct_tau3();

		double h=1e-4;
		CuOSMatrix G_kx_ph_1 = totalGreen_(kx + h,ky).inv();
		CuOSMatrix G_kx_mh_1 = totalGreen_(kx - h,ky).inv();
		CuOSMatrix diffinvG = (G_kx_ph_1 - G_kx_mh_1)/(2.*h);
		for(int oi = 0;oi<3;oi++){
			for(int oj = 0;oj<3;oj++){
				diffinvG(3+oi,oj) = 0;
				diffinvG(oi,3+oj) = 0;
			}
		}
		//Integrated form
		CuOSMatrix inter;
		CuOSMatrix temp = first_d_x*green;			inter = temp*diffinvG;			temp = inter*green;
		CuOSMatrix temp2 = tau3*first_d_x; 	inter = temp2*green; temp2 = inter*tau3;	inter = temp2*diffinvG;	temp2 = inter*green;
		CuOSMatrix result = temp - temp2;

		return result.trace();
	};
private:
	double const tpd_;
	double const tpp_;
	double const tppp_;
	double const ep_;
	std::complex<double> z_;
	RCuOLatticeGreenPeriodized totalGreen_;

	CuOSMatrix construct_tau3() {
		CuOSMatrix tau3 = 0;
		for(int oi = 0; oi < 3; ++oi){
	        tau3(oi,oi) = 1;
			tau3(oi + 3,oi + 3) = -1;
		}
		return tau3;
	}
	CuOSMatrix construct_first_d_x(double kx, double ky) {
		CuOSMatrix result = 0;
		double coskx = std::cos(kx);
		double sinkx = std::sin(kx);
		std::complex<double> exp_kx(coskx, sinkx);
		std::complex<double> exp_mkx = std::conj(exp_kx);
		
		double cosky = std::cos(ky);
		double sinky = std::sin(ky);
		std::complex<double> exp_ky(cosky, sinky);
		std::complex<double> exp_mky = std::conj(exp_ky);
		std::complex<double> I(0.0,1.0);

		//Normal fermi velocity part
		CuOMatrix Vkx=CuOMatrix::Diag(0);Vkx(0, 1) += I*tpd_*exp_mkx;				Vkx(0, 2) += 0;
		Vkx(1, 0) += -I*tpd_*exp_kx;	Vkx(1, 1) += -2.*tppp_*sinkx;				Vkx(1, 2) += -I*tpp_*exp_kx*(1. - exp_mky);
		Vkx(2, 0) += 0;  				Vkx(2, 1) += I*tpp_*exp_mkx*(1. - exp_ky);	Vkx(2, 2) += 0;

		CuOMatrix Vmkx=CuOMatrix::Diag(0);Vmkx(0, 1) += I*tpd_*exp_kx;				Vmkx(0, 2) += 0;
		Vmkx(1, 0) += -I*tpd_*exp_mkx;	Vmkx(1, 1) += 2.*tppp_*sinkx;				Vmkx(1, 2) += -I*tpp_*exp_mkx*(1. - exp_ky);
		Vmkx(2, 0) += 0;  				Vmkx(2, 1) += I*tpp_*exp_kx*(1. - exp_mky);	Vmkx(2, 2) += 0;

		/*
		Vkx(0,0) +=0;						Vkx(0, 1) += 1./2.*I*tpd_*(1. - exp_mkx);				Vkx(0, 2) += 0;
		Vkx(1, 0) += -1./2.*I*tpd_*(1.-exp_kx);	Vkx(1, 1) += 0;											Vkx(1, 2) += - 1./2.*I*tpp_*(1. - exp_kx)*(1. - exp_mky);
		Vkx(2, 0) += 0;  						Vkx(2, 1) += 1./2.*I*tpp_*(1. - exp_mkx)*(1. - exp_ky);	Vkx(2, 2) += 0;

		Vmkx(0,0)+=0;							Vmkx(0, 1) += 1./2.*I*tpd_*(1. - exp_kx);				Vmkx(0, 2) += 0;
		Vmkx(1, 0) += -1./2.*I*tpd_*(1. - exp_mkx);	Vmkx(1, 1) += 0;										Vmkx(1, 2) += - 1./2.*I*tpp_*(1. - exp_mkx)*(1. - exp_ky);
		Vmkx(2, 0) += 0;  							Vmkx(2, 1) += 1./2.*I*tpp_*(1. - exp_kx)*(1. - exp_mky);Vmkx(2, 2) += 0;
		 */
		for(int oi = 0; oi < 3; ++oi){
			for(int oj = 0; oj < 3; ++oj) {
				result(oi,oj) += Vkx(oi, oj);
				result(oi + 3,oj + 3) += std::conj(Vmkx(oi, oj));
			}
		}
		return result;
	}
};


#endif





