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
		double const EkinFM = 2*(ep - 2.*tpp);                                                               // px + py
		double const EkinSM = 4.*tpd*tpd                                                                     // d 
		                    + 2.*((ep - 2.*tpp)*(ep - 2.*tpp) + 2.*tpd*tpd + 6.*tpp*tpp - mu*(ep - 2.*tpp)); //px + py
		 
		double ekin = EkinFM/2. - EkinSM*beta/4.; 
		
		std::ifstream selfFile((dataFolder + "self" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pxgreenFile((dataFolder + "pxgreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pygreenFile((dataFolder + "pygreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pxygreenFile((dataFolder + "pxygreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream pfilling(dataFolder + "pn.dat", std::ios::out | std::ios::app);
		std::ofstream ekinFile(dataFolder + "ekin.dat", std::ios::out | std::ios::app);
		
		std::ofstream greenFile((dataFolder + "dgreen" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		
		//std::ofstream diffFile("diff.dat");
		
		Int::EulerMaclaurin2D<RCuOMatrix> integrator(1.e-10, 4, 12);
		
		for(std::size_t n = 0; n < NMat; ++n) {
			
            std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
			
			double s00_R, s00_I, s01_R, s01_I, s11_R, s11_I;
			selfFile >> dummy >> s00_R >> s00_I >> s01_R >> s01_I >> s11_R >> s11_I;
			std::complex<double> s00(s00_R, s00_I);
			std::complex<double> s01(s01_R, s01_I);
			std::complex<double> s11(s11_R, s11_I);



			RCuMatrix selfEnergy;
			selfEnergy("d_0", "d_0") = s00; selfEnergy("d_0", "d_1") = s01; selfEnergy("d_0", "d_2") = s11; selfEnergy("d_0", "d_3") = s01;
			selfEnergy("d_1", "d_0") = s01; selfEnergy("d_1", "d_1") = s00; selfEnergy("d_1", "d_2") = s01; selfEnergy("d_1", "d_3") = s11;
			selfEnergy("d_2", "d_0") = s11; selfEnergy("d_2", "d_1") = s01; selfEnergy("d_2", "d_2") = s00; selfEnergy("d_2", "d_3") = s01;
			selfEnergy("d_3", "d_0") = s01; selfEnergy("d_3", "d_1") = s11; selfEnergy("d_3", "d_2") = s01; selfEnergy("d_3", "d_3") = s00;
			
			/* Now we do some computation */
			RCuOLatticeGreen latticeGreenRCuO(iomega + mu, tpd, tpp,tppp, ep, selfEnergy);			
			RCuOMatrix green = integrator(latticeGreenRCuO, M_PI/2., M_PI/2.);
			/* No more computation */
			pxgreenFile << iomega.imag() << " "  
			            << green("px_0", "px_0").real() << " " << green("px_0", "px_0").imag()<< " " 
			            << green("px_0", "px_1").real() << " " << green("px_0", "px_1").imag() << " "
			            << green("px_0", "px_2").real() << " " << green("px_0", "px_2").imag() << " "
			            << green("px_0", "px_3").real() << " " << green("px_0", "px_3").imag() << std::endl;
			
			pygreenFile << iomega.imag() << " "  
			            << green("py_0", "py_0").real() << " " << green("py_0", "py_0").imag()<< " " 
			            << green("py_0", "py_1").real() << " " << green("py_0", "py_1").imag() << " "
			            << green("py_0", "py_2").real() << " " << green("py_0", "py_2").imag() << " "
			            << green("py_0", "py_3").real() << " " << green("py_0", "py_3").imag() << std::endl;
			
			pxygreenFile << iomega.imag() << " "  
			             << green("px_0", "py_0").real() << " " << green("px_0", "py_0").imag()<< " " 
			             << green("px_0", "py_1").real() << " " << green("px_0", "py_1").imag() << " " 
			             << green("px_0", "py_2").real() << " " << green("px_0", "py_2").imag() << " " 
						 << green("px_0", "py_3").real() << " " << green("px_0", "py_3").imag() << std::endl;
			
			greenFile << iomega.imag() << " "  
			             << green("d_0", "d_0").real() << " " << green("d_0", "d_0").imag() << " " 
			             << green("d_0", "d_1").real() << " " << green("d_0", "d_1").imag() << " " 
			             << green("d_0", "d_2").real() << " " << green("d_0", "d_2").imag() << " "
			             << green("d_0", "d_3").real() << " " << green("d_0", "d_3").imag() << std::endl;
			
			double temp = .0;
			
			temp += green("px_0", "px_0").real();
			temp += green("px_1", "px_1").real();
			temp += green("px_2", "px_2").real();
			temp += green("px_3", "px_3").real();
			temp += green("py_0", "py_0").real();
			temp += green("py_1", "py_1").real();
			temp += green("py_2", "py_2").real();
			temp += green("py_3", "py_3").real();
			
			np += 2./beta*(temp/8. - (xp/(iomega - xp) - xm/(iomega - xm))/D).real();

			/*here some more computation*/
			RCuOLatticeKineticEnergy latticeKineticEnergyRCuO(iomega + mu, tpd, tpp, tppp,ep, selfEnergy);			
			
			ekin += 2./beta*(integrator(latticeKineticEnergyRCuO, M_PI/2., M_PI/2.).trace()/4. - EkinFM/iomega - EkinSM/(iomega*iomega)).real();
			/* And it's over */
			//std::complex<double> diff = integrator(latticeKineticEnergyRCuO, M_PI/2., M_PI/2.).trace()/4. - EkinFM/iomega - EkinSM/(iomega*iomega); 			
			//diffFile << iomega.imag() << " " << diff.real() << " " << diff.imag() << std::endl;
		}
		
		
		
		if(!selfFile) 
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
		
		//diffFile.close();
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