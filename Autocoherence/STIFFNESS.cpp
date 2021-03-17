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
		if(argc != 6) 
			throw std::runtime_error("Usage : GFULL inputFolder outputFolder dataFolder filename iteration");
		
		std::string inputFolder = argv[1]; 
		std::string outputFolder = argv[2];
		std::string dataFolder = argv[3];
		std::string name = argv[4];
		int const iteration = std::atoi(argv[5]);

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
		std::ofstream stiff(dataFolder + "stiffness.dat", std::ios::out | std::ios::app);
		//
		std::complex<double> stiffness = 0;
		std::complex<double> last_stiffness = 0;
		double error = 1e-2;
		double min_value = 1e-4;
		Int::EulerMaclaurin2D<std::complex<double>> integrator(1.e-2, 4, 12);
		std::size_t n_max = 0;
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

			

			SuperfluidStiffness superfluid(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);
			last_stiffness = integrator(superfluid, M_PI, M_PI);
			stiffness += last_stiffness;
			n_max = n;
			if( std::abs(last_stiffness)/std::abs(stiffness) * (NMat - n) < error/10 || std::abs(stiffness) < min_value){
				break;
			}
		}
	//Because we sumed on the Matsubara Frequencies, we must normalize by a $\beta$ factor (this is in the formula of course)
	stiffness/=beta;
	std::cout << "This operation went to the " << n_max << "th Matsubara frequency" << std::endl;
    	stiff << iteration << " " << stiffness.real() << " " << stiffness.imag() << std::endl;
    	stiff.close();
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












