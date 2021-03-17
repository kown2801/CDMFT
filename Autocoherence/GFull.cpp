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
		if(argc != 6) 
			throw std::runtime_error("Usage : GFULL inputFolder outputFolder dataFolder filename iteration");
		
		std::string inputFolder = argv[1]; 
		std::string outputFolder = argv[2];
		std::string dataFolder = argv[3];
		std::string name = argv[4];
		int const iteration = std::atoi(argv[5]);

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
		
                
        /***************************/
        /* We read the Link file   */
        json_spirit::mArray jLink;       
		std::size_t nSite_;
        newIO::readLinkFromParams(jLink, nSite_, outputFolder, readParams);
        /***************************/
        /*****************************************************************/
        /* Now we create the map object that will contain the components */
        std::map<std::string,std::complex<double> > component_map;
        std::map<std::string,std::vector<std::pair<std::size_t,std::size_t> > > inverse_component_map;
        for(std::size_t i=0;i<jLink.size();i++){
            for(std::size_t j=0;j<jLink.size();j++){
                component_map[jLink[i].get_array()[j].get_str()] = 0;
                if ( inverse_component_map.find(jLink[i].get_array()[j].get_str()) == inverse_component_map.end() ) {
                    inverse_component_map[jLink[i].get_array()[j].get_str()] = std::vector<std::pair<std::size_t,std::size_t> >();
                }
                inverse_component_map[jLink[i].get_array()[j].get_str()].push_back(std::pair<std::size_t,std::size_t>(i,j));
            }
        }
        /*****************************************************************/
        newIO::GenericReadFunc readSelf(dataFolder + "self" + std::to_string(iteration) + ".json","");

		if(!readSelf.good()){ 
			throw std::runtime_error("Error while reading selfenergy file.");
		}

		std::ofstream pxgreenFile((dataFolder + "pxgreen" + std::to_string(iteration) + ".dat").c_str());
		std::ofstream pygreenFile((dataFolder + "pygreen" + std::to_string(iteration) + ".dat").c_str());
		std::ofstream pxygreenFile((dataFolder + "pxygreen" + std::to_string(iteration) + ".dat").c_str());
		std::ofstream pfilling(dataFolder + "pn.dat", std::ios::out | std::ios::app);
		std::ofstream ekinFile(dataFolder + "ekin.dat", std::ios::out | std::ios::app);	
		std::ofstream greenFile((dataFolder + "dgreen" + std::to_string(iteration) + ".dat").c_str());
		//
		
		Int::EulerMaclaurin2D<RCuOMatrix> integrator(1.e-4, 4, 12);
		
		for(std::size_t n = 0; n < NMat; ++n) {
			
			std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
			/***********************************************/
			/* We load the different selfEnergy components */
 			for (auto &p : component_map)
            {
                if(p.first != "empty"){ //We don't read the empty component
                    p.second = readSelf(p.first)->getFunction(n);
                }
            } 
			/***********************************************/
			RCuMatrix selfEnergy;
			self_constraints(component_map);
			newIO::component_map_to_matrix(jLink,selfEnergy,component_map);
			
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
			//The down spins contribute with a minus sign because of the Nambu convention
			temp -= green("px_0Down", "px_0Down").real();
			temp -= green("px_1Down", "px_1Down").real();
			temp -= green("px_2Down", "px_2Down").real();
			temp -= green("px_3Down", "px_3Down").real();
			temp -= green("py_0Down", "py_0Down").real();
			temp -= green("py_1Down", "py_1Down").real();
			temp -= green("py_2Down", "py_2Down").real();
			temp -= green("py_3Down", "py_3Down").real();
			
			//np is the number of electrons per site per spin

			np += 2./beta*(temp/16 - (xp/(iomega - xp) - xm/(iomega - xm))/D).real();
			
			RCuOLatticeKineticEnergy latticeKineticEnergyRCuO(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);			
			
			ekin += 2./beta*(integrator(latticeKineticEnergyRCuO, M_PI/2., M_PI/2.).trace()/8. - EkinFM/iomega - EkinSM/(iomega*iomega)).real();		
		}
		
		
		pfilling << iteration << " " << np << std::endl;
		ekinFile << iteration << " " << ekin << std::endl;
		
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












