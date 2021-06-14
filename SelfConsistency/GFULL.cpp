#include "Patrick/Integrators.h"
#include "Patrick/Hyb.h"
#include "Patrick/Plaquette/Plaquette.h"
#include "IO.h"

double fermi(double arg) {
	return arg > .0 ? std::exp(-arg)/(1. + std::exp(-arg)) : 1./(1. + std::exp(arg));
};
void addComponentToObject(const RCuOMatrix& matrix, const std::string left, const std::string right, json& jObject,double const beta){
	std::string index = left + "," + right;
	if(!jObject[index].count("real")){
		jObject[index]["real"] = json::array();
		jObject[index]["imag"] = json::array();
		jObject[index]["dataType"] = "FermionicMatsubaraFrequencies";
		jObject[index]["beta"] = beta;
	}
	jObject[index]["real"].push_back(matrix(left,right).real());
	jObject[index]["imag"].push_back(matrix(left,right).imag());
}
int main(int argc, char** argv)
{
	try {
		if(argc != 6) 
			throw std::runtime_error("Usage : GFULL inputFolder outputFolder dataFolder inputFilename iteration");
		
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
		IO::readJsonFileAtIteration(inputFolder + name,".meas.json",iteration,jParams);
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
        IO::readJsonFileAtIteration(dataFolder + "self",".json",iteration,jSelf);
        Hyb::accountForDifferentBeta(jSelf,beta);
		
		json jPx,jPy,jPxy,jD;
		
		Int::EulerMaclaurin2D<RCuOMatrix> integrator(1.e-4, 4, 12);
		
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
			
			RCuOLatticeGreen latticeGreen(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);			
			RCuOMatrix green = integrator(latticeGreen, M_PI/2., M_PI/2.);
			/********************************************************************************/
			/* We save some of the components of the Green's function to file for later use */
			addComponentToObject(green,"px_0Up", "px_0Up",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_1Up",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_2Up",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_3Up",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_0Down",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_1Down",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_2Down",jPx,beta);
			addComponentToObject(green,"px_0Up", "px_3Down",jPx,beta);
			
			addComponentToObject(green,"py_0Up", "py_0Up",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_1Up",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_2Up",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_3Up",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_0Down",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_1Down",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_2Down",jPy,beta);
			addComponentToObject(green,"py_0Up", "py_3Down",jPy,beta);
			
			addComponentToObject(green,"px_0Up", "py_0Up",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_1Up",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_2Up",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_3Up",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_0Down",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_1Down",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_2Down",jPxy,beta);
			addComponentToObject(green,"px_0Up", "py_3Down",jPxy,beta);
			
			addComponentToObject(green,"d_0Up", "d_0Up",jD,beta);
			addComponentToObject(green,"d_0Up", "d_1Up",jD,beta);
			addComponentToObject(green,"d_0Up", "d_2Up",jD,beta);
			addComponentToObject(green,"d_0Up", "d_1Down",jD,beta);
			addComponentToObject(green,"d_1Up", "d_2Down",jD,beta);
			/********************************************************************************/
			
			/************************************/
			/* We compute the oxygen occupation */
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
			
			np += 2./beta*(temp/16 - (xp/(iomega - xp) - xm/(iomega - xm))/D).real();
			//np is the number of electrons per site per spin
			/************************************/
			
			/*********************************/
			/* We compute the kinetic energy */
			RCuOLatticeKineticEnergy latticeKineticEnergyRCuO(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);			
			
			ekin += 2./beta*(integrator(latticeKineticEnergyRCuO, M_PI/2., M_PI/2.).trace()/8. - EkinFM/iomega - EkinSM/(iomega*iomega)).real();		
			/*********************************/
		}
		IO::writeInJsonDataFile(dataFolder + "pxgreen.json",iteration,jPx);
		IO::writeInJsonDataFile(dataFolder + "pygreen.json",iteration,jPy);
		IO::writeInJsonDataFile(dataFolder + "pxygreen.json",iteration,jPxy);
		IO::writeInJsonDataFile(dataFolder + "dgreen.json",iteration,jD);
		
		std::ofstream pfilling(dataFolder + "np.dat", std::ios::out | std::ios::app);
		std::ofstream ekinFile(dataFolder + "ekin.dat", std::ios::out | std::ios::app);	
		
		pfilling << iteration << " " << np << std::endl;
		ekinFile << iteration << " " << ekin << std::endl;
		
		pfilling.close();
		ekinFile.close();
		
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












