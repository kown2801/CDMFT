#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"

//Vielleicht noch dritter moment der gp, aber wie am besten .... ?

double fermi(double arg) {
	return arg > .0 ? std::exp(-arg)/(1. + std::exp(-arg)) : 1./(1. + std::exp(arg));
};

int main(int argc, char** argv)
{
	try {
		if(argc != 5) 
			throw std::runtime_error("Usage : GFULL inputFolder dataFolder filename iteration");
		
		std::string inputFolder = argv[1]; 
		std::string dataFolder = argv[2];
		std::string name = argv[3];
		int const iteration = std::atoi(argv[4]);

		newIO::GenericReadFunc readParams(inputFolder + name + boost::lexical_cast<std::string>(iteration) + ".meas.json","Parameters");

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
		
		double const A = ep - 2.*tpp - mu;
		double const B = 2.*tpd*tpd + 4.*tpp*tpp + 2.*tppp*tppp;		
		double const D = std::sqrt(A*A + 4.*B);
		
		double const xp = (A + D)/2.;
		double const xm = (A - D)/2.;
		
		std::string dummy;
		double np = (xp*fermi(beta*xp) - xm*fermi(beta*xm))/D;
		
		//Per unit cell
		double const EkinFM = 2*(ep - 2.*tpp);                                                               			//px + py
		double const EkinSM = 4.*tpd*tpd                                                                     				//d 
						    + 2.*((ep - 2.*tpp)*(ep - 2.*tpp) + 2.*tpd*tpd + 4.*tpp*tpp + 2.*tppp*tppp - mu*(ep - 2.*tpp)); //px + py
		
		double ekin = EkinFM/2. - EkinSM*beta/4.;
		
		std::ifstream selfFile((dataFolder + "self" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pxgreenFile((dataFolder + "pxgreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pygreenFile((dataFolder + "pygreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pxygreenFile((dataFolder + "pxygreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pfilling(dataFolder + "pn.dat", std::ios::out | std::ios::app);
		std::ofstream ekinFile(dataFolder + "ekin.dat", std::ios::out | std::ios::app);	
		std::ofstream greenFile((dataFolder + "dgreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		//
		
		Int::EulerMaclaurin2D<RCuOMatrix> integrator(1.e-4, 4, 12);
		
		for(std::size_t n = 0; n < NMat; ++n) {
			
			std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
						
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
			
			RCuOLatticeGreen latticeGreen(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);			
			RCuOMatrix green = integrator(latticeGreen, M_PI/2., M_PI/2.);
			
			pxgreenFile << iomega.imag() << " "  
			            << green("px_0Up", "px_0Up").real() << " " << green("px_0Up", "px_0Up").imag() << " " 
			            << green("px_0Up", "px_1Up").real() << " " << green("px_0Up", "px_1Up").imag() << " " 
			            << green("px_0Up", "px_2Up").real() << " " << green("px_0Up", "px_2Up").imag() << " " 
			            << green("px_0Up", "px_3Up").real() << " " << green("px_0Up", "px_3Up").imag() << " "
			            << green("px_0Up", "px_0Down").real() << " " << green("px_0Up", "px_0Down").imag() << " " 
			            << green("px_0Up", "px_1Down").real() << " " << green("px_0Up", "px_1Down").imag() << " "
			            << green("px_0Up", "px_2Down").real() << " " << green("px_0Up", "px_2Down").imag() << " "
			            << green("px_0Up", "px_3Down").real() << " " << green("px_0Up", "px_3Down").imag() << std::endl;
			
			pygreenFile << iomega.imag() << " "  
						<< green("py_0Up", "py_0Up").real() << " " << green("py_0Up", "py_0Up").imag() << " " 
			            << green("py_0Up", "py_1Up").real() << " " << green("py_0Up", "py_1Up").imag() << " " 
			            << green("py_0Up", "py_2Up").real() << " " << green("py_0Up", "py_2Up").imag() << " " 
						<< green("py_0Up", "py_3Up").real() << " " << green("py_0Up", "py_3Up").imag() << " "
			            << green("py_0Up", "py_0Down").real() << " " << green("py_0Up", "py_0Down").imag() << " " 
			            << green("py_0Up", "py_1Down").real() << " " << green("py_0Up", "py_1Down").imag() << " "
						<< green("py_0Up", "py_2Down").real() << " " << green("py_0Up", "py_2Down").imag() << " "
						<< green("py_0Up", "py_3Down").real() << " " << green("py_0Up", "py_3Down").imag() << std::endl;
			
			pxygreenFile << iomega.imag() << " "  
			             << green("px_0Up", "py_0Up").real() << " " << green("px_0Up", "py_0Up").imag() << " " 
			             << green("px_0Up", "py_1Up").real() << " " << green("px_0Up", "py_1Up").imag() << " " 
			             << green("px_0Up", "py_2Up").real() << " " << green("px_0Up", "py_2Up").imag() << " " 
			             << green("px_0Up", "py_3Up").real() << " " << green("px_0Up", "py_3Up").imag() << " " 
			             << green("px_0Up", "py_0Down").real() << " " << green("px_0Up", "py_0Down").imag() << " " 
						 << green("px_0Up", "py_1Down").real() << " " << green("px_0Up", "py_1Down").imag() << " "
			             << green("px_0Up", "py_2Down").real() << " " << green("px_0Up", "py_2Down").imag() << " "
			             << green("px_0Up", "py_3Down").real() << " " << green("px_0Up", "py_3Down").imag() << std::endl;
			
			
			greenFile << iomega.imag() << " "  
			          << green("d_0Up", "d_0Up").real() << " " << green("d_0Up", "d_0Up").imag() << " " 
			          << green("d_0Up", "d_1Up").real() << " " << green("d_0Up", "d_1Up").imag() << " " 
			          << green("d_0Up", "d_2Up").real() << " " << green("d_0Up", "d_2Up").imag() << " " 
			          << green("d_0Up", "d_1Down").real() << " " << green("d_0Up", "d_1Down").imag() << " "
			          << green("d_1Up", "d_2Down").real() << " " << green("d_1Up", "d_2Down").imag() << std::endl;
			
			
			double temp = .0;
			
			temp += green("px_0Up", "px_0Up").real();
			temp += green("px_1Up", "px_1Up").real();
			temp += green("px_2Up", "px_2Up").real();
			temp += green("px_3Up", "px_3Up").real();
			temp += green("py_0Up", "py_0Up").real();
			temp += green("py_1Up", "py_1Up").real();
			temp += green("py_2Up", "py_2Up").real();
			temp += green("py_3Up", "py_3Up").real();
			
			np += 2./beta*(temp/8. - (xp/(iomega - xp) - xm/(iomega - xm))/D).real();
			
			RCuOLatticeKineticEnergy latticeKineticEnergyRCuO(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);			
			
			ekin += 2./beta*(integrator(latticeKineticEnergyRCuO, M_PI/2., M_PI/2.).trace()/8. - EkinFM/iomega - EkinSM/(iomega*iomega)).real();			
		}
		
		if(selfFile >> dummy) 
			throw std::runtime_error("Error while reading selfenergy file.");
		
		pfilling << iteration << " " << np << std::endl;
		ekinFile << iteration << " " << ekin << std::endl;
		
		selfFile.close();
		pxgreenFile.close();
		pygreenFile.close();
		pxygreenFile.close();
		pfilling.close();
		ekinFile.close();
		
		greenFile.close();
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












