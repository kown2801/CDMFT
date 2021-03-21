#include "Patrick/Integrators.h"
#include "Patrick/Hyb.h"
#include "Patrick/Plaquette/Plaquette.h"
#include "IO.h"

double fermi(double arg) {
	return arg > .0 ? std::exp(-arg)/(1. + std::exp(-arg)) : 1./(1. + std::exp(arg));
};
int main(int argc, char** argv)
{
	try {
		if(argc != 6) 
			throw std::runtime_error("Usage : GFULL inputFolder outputFolder dataFolder filename iteration");
		
		/******************************/
        /* Initialisation of variable */
        std::string inputFolder = argv[1];
        std::string outputFolder = argv[2];
        std::string dataFolder = argv[3];
        std::string name = argv[4];
        int const iteration = std::atoi(argv[5]);
		/******************************/
        /*   Reading the parameters   */
        /******************************/
        json jParams;
		IO::readJsonFile(inputFolder + name + std::to_string(iteration) + ".meas.json",jParams);
		jParams = jParams["Parameters"];

		double const mu = jParams["mu"];
        double const beta = jParams["beta"];
        double const tpd = jParams["tpd"];
        double const tpp = jParams["tpp"];
        double const ep = jParams["ep"];
        double tppp = tpp;
		if(exists(jParams,"tppp")){
			tppp = jParams["tppp"];
            std::cout << "We have tppp different than tpp" << std::endl;
		}
		double const EGreen = jParams["EGreen"];
		unsigned int const NMat = beta*EGreen/(2*M_PI) + 1;
	
		/***************************/
        /* We read the Link file   */
        json jLink;       
        IO::readLinkFromParams(jLink, outputFolder, jParams);
        /***************************/
       /*****************************************************************/
        /* Now we create the map object that will contain the components */
        std::map<std::string,std::complex<double> > component_map;
        std::map<std::string,std::vector<std::pair<std::size_t,std::size_t> > > inverse_component_map;
        for(std::size_t i=0;i<jLink.size();i++){
            for(std::size_t j=0;j<jLink.size();j++){
                component_map[jLink[i][j]] = 0;
                if ( inverse_component_map.find(jLink[i][j]) == inverse_component_map.end() ) {
                    inverse_component_map[jLink[i][j]] = std::vector<std::pair<std::size_t,std::size_t> >();
                }
                inverse_component_map[jLink[i][j]].push_back(std::pair<std::size_t,std::size_t>(i,j));
            }
        }
        /*****************************************************************/

        json jSelf;
        IO::readJsonFile(dataFolder + "self" + std::to_string(iteration) + ".json",jSelf);
        Hyb::accountForDifferentBeta(jSelf,beta);

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
                    p.second = std::complex<double>(jSelf[p.first]["real"][n],jSelf[p.first]["imag"][n]);
                }
            } 
			/***********************************************/
			RCuMatrix selfEnergy;
			self_constraints(component_map);
			IO::component_map_to_matrix(jLink,selfEnergy,component_map);

			

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












