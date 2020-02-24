#ifndef __UTILITIES
#define __UTILITIES

#include <vector>
#include <complex>
#include <boost/random.hpp>

//Input fŸr Matrizen ist Row-Major
//BOOST_STATIC_ASSERT(numeric_limits<double>::is_iec559) for memset ?????

extern "C" {
	double dasum_(int const*, double const*, int const*);
	double dnrm2_(int const*, double const*, int const*);
	double ddot_(const int*, const double*, const int*, const double*, const int*);
	void   dswap_(const int*, double*, int const*, double*, int const*);
	void   daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
	void   dger_(const int*, const int*, double const*, double const*, int const*, double const*, int const*, double*, int const*);
	void   dscal_(int const*, double const*, double*, int const*);
	void   dgemm_(const char*, const char*, int const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
	void   dcopy_(int const*, double const*, int const* , double*, int const*);
	void   dgemv_(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void   dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);
	void   zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
}

namespace Ut {		
	//-----------------------------------------------------------------------------------------------------------------------------------------------------
	
	typedef std::complex<double> complex;
	typedef boost::mt19937 EngineType;
	typedef boost::uniform_real<double> UniformDistribution;
	typedef boost::variate_generator<EngineType&, UniformDistribution> UniformRngType;
	
	//-------------------------------------------------------------------------------------------------------------------------------------------------------
	
	template<class T>
	struct MeasEntry {
		MeasEntry() {};
		MeasEntry(T mean, T error) : mean_(mean), error_(error) {};
		T mean() const { return mean_;};
		T error() const { return error_;};
	private:
		T mean_; T error_;
	};
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------
	
	template<class T>
	void read(std::ifstream& file, T& t) {
		file.read(reinterpret_cast<char*>(&t), sizeof(T));
	};
	
	template<class T>
	void write(std::ofstream& file, T t) {
		file.write(reinterpret_cast<char*>(&t), sizeof(T));
	};
};


#endif