#ifndef __INTEGRATORS
#define __INTEGRATORS

#include <cmath>
#include <iostream>
#include <valarray>
#include <vector>
#include <utility>
#include <set>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <cmath>
#include <limits>


namespace Int {
	template<class RET>
	class Gauss2D {
	public:
		Gauss2D(double error, int dephtMin, int dephtMax) : 
		error_(error), 
		dephtMin_(dephtMin), 
		dephtMax_(dephtMax),
		v_(10./49.),
		w_(9./196.),
		h_(std::sqrt(7./15.)),
		k_(std::sqrt(7.)/3.) {
			if(dephtMin < 1) throw std::runtime_error("Gauss2D: iterMin > 0 !");
			if(dephtMin_ > dephtMax_) throw std::runtime_error("Gauss2D: iterMin <= iterMax !"); 
		};
		
		template<typename Func>
		RET const& operator()(Func& func, double X, double Y) {
			RET workSpace[4*(dephtMax_ + 1)]; 		
			result_ = .0;  Rec(result_, func, .0, .0, X, Y, std::numeric_limits<double>::max(), dephtMin_, workSpace);
			result_ = .0;  Rec(result_, func, .0, .0, X, Y, error_*result_.norm(), dephtMax_, workSpace);
			return result_;
		};
	private:
		double const error_;
		int const dephtMin_;
		int const dephtMax_;
		double const v_;
		double const w_;
		double const h_;
		double const k_;
		
		RET temp_;
		RET result_;		
		
		template<typename Func>
		void Rec(RET& ret, Func& func, double x, double y, double X, double Y, double error, int countDown, RET* workSpace) {
			Cub(workSpace[0], func, x + X/2., y + Y/2., X/2., Y/2.);
			Cub(workSpace[1], func, x + X/2., y - Y/2., X/2., Y/2.);
			Cub(workSpace[2], func, x - X/2., y + Y/2., X/2., Y/2.);
			Cub(workSpace[3], func, x - X/2., y - Y/2., X/2., Y/2.);
			
			temp_ += workSpace[0];
			temp_ += workSpace[1];
			temp_ += workSpace[2];
			temp_ += workSpace[3]; 
			temp_ *= .25;
			ret -= temp_;
			if(ret.norm() < error || countDown == 0) 
				ret = temp_;
			else {
				ret = .0;
				Rec(workSpace[0], func, x + X/2., y + Y/2., X/2., Y/2., 2.*error, countDown - 1, workSpace + 4); ret += workSpace[0];
				Rec(workSpace[1], func, x + X/2., y - Y/2., X/2., Y/2., 2.*error, countDown - 1, workSpace + 4); ret += workSpace[1];
				Rec(workSpace[2], func, x - X/2., y + Y/2., X/2., Y/2., 2.*error, countDown - 1, workSpace + 4); ret += workSpace[2];
				Rec(workSpace[3], func, x - X/2., y - Y/2., X/2., Y/2., 2.*error, countDown - 1, workSpace + 4); ret += workSpace[3];
				ret *= .25;
			}
		}
		
		template<typename Func>
		void Cub(RET& ret, Func& func, double x, double y, double X, double Y) {
			ret = func(x + h_*X, y);
			ret += func(x - h_*X, y);
			ret += func(x, y + h_*Y);
			ret += func(x, y - h_*Y);
			ret *= v_;

			temp_ = func(x + k_*X, y + k_*Y);
			temp_ += func(x - k_*X, y + k_*Y);
			temp_ += func(x + k_*X, y - k_*Y);
			temp_ += func(x - k_*X, y - k_*Y);

			ret += w_*temp_;
		};
	};
	
	
	
	template<class RET>
	struct EulerMaclaurin2D {
		EulerMaclaurin2D(double error, int nMin, int nMax) : error_(error), nMin_(nMin), nMax_(nMax) {
			if(nMin_ < 1) throw std::runtime_error("EulerMaclaurin2D: nMin > 0 !");
			if(nMin_ > nMax_) throw std::runtime_error("EulerMaclaurin2D: nMin <= nMax !"); 
		};
		template<typename Func>
		RET const& operator()(Func& func, double X, double Y) {
			int steps = 1 << (nMin_ - 1);
			RET error; result_ = .0;
			do {
				steps *= 2;				
				error = result_; 
				
				result_ = .0;
				for(int i = 0; i < steps; i++) 
					for(int j = 0; j < steps; j++) {
						double kx = -X*(2.*i - steps + 1.)/static_cast<double>(steps - 1);
						double ky = -Y*(2.*j - steps + 1.)/static_cast<double>(steps - 1);
						double fact(1.);
						fact *= i%(steps - 1) == 0 ? .5 : 1.;
						fact *= j%(steps - 1) == 0 ? .5 : 1.;
						
						result_ += fact*func(kx, ky);
					}
				
				result_ *= 1./((steps - 2.)*(steps - 2.) + 2.*(steps - 2.) + 1.);
				error -= result_;
		    } while(abs(error) > error_*abs(result_) && steps < (1 << nMax_));
			
			return result_;
		};
	private:
		double const error_;
		int const nMin_;
		int const nMax_;
		
		RET result_;
	};
};

#endif