#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"

double fermi(double arg) {
	return arg > .0 ? std::exp(-arg)/(1. + std::exp(-arg)) : 1./(1. + std::exp(arg));
};
int main(int argc, char** argv)
{
	try {
		if(argc != 5) 
			throw std::runtime_error("Usage : SUSCEPTIBILITY inputFolder dataFolder filename iteration");
 		std::string inputFolder = argv[1]; 
		std::string dataFolder = argv[2];
		std::string name = argv[3];
		int const iteration = std::atoi(argv[4]);

		std::cout << "Iteration : " << iteration << std::endl;
		newIO::GenericReadFunc readParams(inputFolder + name + std::to_string(iteration) + ".meas.json","Parameters");

		double const mu = readParams("mu")->getDouble();
		double const beta = readParams("beta")->getDouble();
		double const tpd = readParams("tpd")->getDouble();
		double const tpp = readParams("tpp")->getDouble();
		double tppp = tpp;
		//We want to read tppp if it is defined in the parameter file
		bool existstppp;
		const newIO::GenericReader* tpppRead = readParams("tppp",existstppp);
		if(tpppRead) {
			tppp = tpppRead->getDouble();
		}
		//End of tppp read
		double const ep = readParams("ep")->getDouble();
		unsigned int const NMat = beta*readParams("EGreen")->getInt()/(2*M_PI) + 1;
	
		std::string dummy;
		std::ifstream selfFile((dataFolder + "self" + std::to_string(iteration) + ".dat").c_str());
		if(!selfFile.is_open()){
			std::cout << "No self-Energy file found. Continuing as if the self-Energy was zero" << std::endl;
		}
		std::ofstream susceptibilityFile(dataFolder + "susceptibility.dat", std::ios::out | std::ios::app);
		//
		CuOMatrix susceptibility0_0 = 0;
		CuOMatrix susceptibilityPI_0 = 0;
		CuOMatrix susceptibilityPI_PI = 0;
		double error = 1e-2;
		double min_value = 1e-4;
		Int::EulerMaclaurin2D<CuOMatrix> integrator(1.e-2, 4, 12);
		std::size_t n_max = 0;
		for(std::size_t n = 0; n < NMat; ++n) {
			
			std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
			//We intialize the self energy						
			double s00_R, s00_I, s01_R, s01_I, s11_R, s11_I,pphi_R, pphi_I, mphi_R, mphi_I;;
			selfFile >> dummy >> s00_R >> s00_I >> s01_R >> s01_I >> s11_R >> s11_I >> pphi_R >> pphi_I >> mphi_R >> mphi_I;
			std::complex<double> s00(s00_R, s00_I);
			std::complex<double> s01(s01_R, s01_I);
			std::complex<double> s11(s11_R, s11_I);
			std::complex<double> pphi(pphi_R, pphi_I);
			std::complex<double> mphi(mphi_R, mphi_I);
			
			RCuMatrix selfEnergy;
			selfEnergy("d_0Up", "d_0Up") = s00; selfEnergy("d_0Up", "d_1Up") = s01; selfEnergy("d_0Up", "d_2Up") = s11; selfEnergy("d_0Up", "d_3Up") = s01;
			selfEnergy("d_1Up", "d_0Up") = s01; selfEnergy("d_1Up", "d_1Up") = s00; selfEnergy("d_1Up", "d_2Up") = s01; selfEnergy("d_1Up", "d_3Up") = s11;
			selfEnergy("d_2Up", "d_0Up") = s11; selfEnergy("d_2Up", "d_1Up") = s01; selfEnergy("d_2Up", "d_2Up") = s00; selfEnergy("d_2Up", "d_3Up") = s01;
			selfEnergy("d_3Up", "d_0Up") = s01; selfEnergy("d_3Up", "d_1Up") = s11; selfEnergy("d_3Up", "d_2Up") = s01; selfEnergy("d_3Up", "d_3Up") = s00;	
			
			selfEnergy("d_0Down", "d_0Down") = -std::conj(s00); selfEnergy("d_0Down", "d_1Down") = -std::conj(s01); selfEnergy("d_0Down", "d_2Down") = -std::conj(s11); selfEnergy("d_0Down", "d_3Down") = -std::conj(s01);
			selfEnergy("d_1Down", "d_0Down") = -std::conj(s01); selfEnergy("d_1Down", "d_1Down") = -std::conj(s00); selfEnergy("d_1Down", "d_2Down") = -std::conj(s01); selfEnergy("d_1Down", "d_3Down") = -std::conj(s11);
			selfEnergy("d_2Down", "d_0Down") = -std::conj(s11); selfEnergy("d_2Down", "d_1Down") = -std::conj(s01); selfEnergy("d_2Down", "d_2Down") = -std::conj(s00); selfEnergy("d_2Down", "d_3Down") = -std::conj(s01);
			selfEnergy("d_3Down", "d_0Down") = -std::conj(s01); selfEnergy("d_3Down", "d_1Down") = -std::conj(s11); selfEnergy("d_3Down", "d_2Down") = -std::conj(s01); selfEnergy("d_3Down", "d_3Down") = -std::conj(s00);	
			
			selfEnergy("d_0Up", "d_1Down") = 
			selfEnergy("d_0Down", "d_1Up") = 
			selfEnergy("d_1Up", "d_0Down") = 
			selfEnergy("d_1Down", "d_0Up") = 
			selfEnergy("d_2Up", "d_3Down") = 
			selfEnergy("d_2Down", "d_3Up") =
			selfEnergy("d_3Up", "d_2Down") = 
			selfEnergy("d_3Down", "d_2Up") = (pphi - mphi).real()/2.; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			selfEnergy("d_1Up", "d_2Down") = 
			selfEnergy("d_1Down", "d_2Up") = 
			selfEnergy("d_2Up", "d_1Down") = 
			selfEnergy("d_2Down", "d_1Up") = 
			selfEnergy("d_3Up", "d_0Down") = 
			selfEnergy("d_3Down", "d_0Up") = 
			selfEnergy("d_0Up", "d_3Down") = 
			selfEnergy("d_0Down", "d_3Up") = (mphi - pphi).real()/2.; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			//We create the susceptibility using the selfEnergy object to be able to integrate on the full Brillouin zone
			
			SpinSusceptibility susceptibilityFunction0_0(iomega + mu, tpd, tpp, tppp, ep, selfEnergy,0,0);
			susceptibility0_0 += integrator(susceptibilityFunction0_0, M_PI, M_PI);
			SpinSusceptibility susceptibilityFunctionPI_0(iomega + mu, tpd, tpp, tppp, ep, selfEnergy,M_PI,0);
			susceptibilityPI_0 += integrator(susceptibilityFunctionPI_0, M_PI, M_PI);
			SpinSusceptibility susceptibilityFunctionPI_PI(iomega + mu, tpd, tpp, tppp, ep, selfEnergy,M_PI,M_PI);
			susceptibilityPI_PI += integrator(susceptibilityFunctionPI_PI, M_PI, M_PI);
			
			susceptibilityFile 	<< iomega.imag() << " " << susceptibility0_0(0,0).real()/beta << " "
    					 	<< susceptibilityPI_0(0,0).real()/beta << " "
    					 	<< susceptibilityPI_PI(0,0).real()/beta << " "
						<< ((susceptibility0_0(1,1) + susceptibility0_0(2,2))/2.).real()/beta << " " 
						<< ((susceptibilityPI_0(1,1) + susceptibilityPI_0(2,2))/2.).real()/beta << " " 
						<< ((susceptibilityPI_PI(1,1) + susceptibilityPI_PI(2,2))/2.).real()/beta << " " 
						<< susceptibilityPI_PI(0,1).real()/beta << std::endl;
			
			/*
			 n_max = n;
			if( std::abs(last_stiffness)/std::abs(stiffness) * (NMat - n) < error/10 || std::abs(stiffness) < min_value){
				break;
			}
			*/
		}
		susceptibility0_0/=beta;
		susceptibilityPI_0/=beta;
		susceptibilityPI_PI/=beta;
	susceptibilityFile.close();
    }
	catch(std::exception& exc) {
		std::cerr << exc.what() << "\n";
		return -1;
	}
	catch(...) {
		std::cerr << "Fatal Error: Unknown Exception!\n";
		return -2;
	}
	return 0;
}












