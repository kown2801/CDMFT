#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"

//Vielleicht noch dritter moment der gp, aber wie am besten .... ?

double fermi(double arg) {
	return arg > .0 ? std::exp(-arg)/(1. + std::exp(-arg)) : 1./(1. + std::exp(arg));
};
void readLinkFile(std::string fileName, json_spirit::mArray& jLink,std::string outputFolder){
	std::ifstream file(outputFolder + fileName); 
	if(file) {
		json_spirit::mValue temp;
		json_spirit::read(file, temp); 
		jLink = temp.get_array();
	}else{
		throw std::runtime_error(outputFolder + fileName + " not found.");
	}
}
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
		//First we need to load the structure of the self file and then relate it to the Link file
		std::ifstream structureFile((dataFolder + "structure" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		structureFile>>dummy; //We don't need the component frequency so we get rid of it
		std::vector<std::string> component_ids;
		std::string component_name;			
	       	while(structureFile >> component_name) {
			component_ids.push_back(component_name);
		}
		json_spirit::mArray jLink;
		readLinkFile(readParams("LINK")->getString(),jLink,outputFolder);
		std::size_t nSite_ = jLink.size()/2;
		for(std::size_t n = 0; n < NMat; ++n) {
			
			std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
						
			//First we need to load all the selfEnergy components into a component_map
			std::map<std::string,std::complex<double> > component_map;
			selfFile >> dummy;//We discard the Matsubara frequency
			for(std::size_t i = 0 ; i < component_ids.size() ; i++){
				double this_component_real,this_component_imag;		
				selfFile >> this_component_real >> this_component_imag;
				component_map[component_ids[i]] = std::complex<double>(this_component_real,this_component_imag);
			}
			RCuMatrix selfEnergy;
			//Now we can initialize the selfEnergy
			for(std::size_t i=0;i<jLink.size();i++){
				for(std::size_t j=0;j<jLink.size();j++){
					if (jLink[i].get_array()[j].get_str() != "empty"){
						std::complex<double> this_component = component_map[jLink[i].get_array()[j].get_str()];		
						selfEnergy(i,j) = this_component;
						//We need to beware to initialize the down component according to the Nambu convention
						if(i >= nSite_ && j>= nSite_){
							selfEnergy(i,j) = -std::conj(selfEnergy(i,j));
						}
					}
				}
			}
			//Then we need to make sure the anormal part is symmtric and without a phase. (apparently this is very important)
			//This is specific to our Link.json file and should be removed or should should the Link File change unfortunately 
			std::complex<double> mphi_component = component_map["mphi"];		
			std::complex<double> pphi_component = component_map["pphi"];		
			for(std::size_t i=0;i<jLink.size();i++){
				for(std::size_t j=0;j<jLink.size();j++){
					if (jLink[i].get_array()[j].get_str() == "pphi"){
						selfEnergy(i,j) = (pphi_component - mphi_component).real()/2;
					}else if (jLink[i].get_array()[j].get_str() == "mphi"){
						selfEnergy(i,j) = (mphi_component - pphi_component).real()/2;
					}
				}
			}
			//Now we compute the Full Cluster Green's function
			RCuOLatticeGreen latticeGreen(iomega + mu, tpd, tpp, tppp, ep, selfEnergy);			
			RCuOMatrix green = integrator(latticeGreen, M_PI/2., M_PI/2.);
			
			//We then output the results we need in the data files
			//I don't really know what we need right now so I don't change this part
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
			
			//We then compute the oxygen occupation (just as in the non-antiferromagnetic case, there is no difference of course)
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












