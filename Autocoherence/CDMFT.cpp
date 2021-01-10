#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"
#include <json_spirit.h>


void readScalSites(std::string obs, newIO::GenericReadFunc& readParams, int iteration,std::string outputFolder) {
	std::ofstream file((outputFolder + obs + "Sites.dat").c_str(), std::ios_base::out | std::ios_base::app);				
	file << iteration; 
	
	for(int i = 0; i < 4; ++i) {
		std::string s = boost::lexical_cast<std::string>(i); 
		file << " " << readParams(obs + "_" + s)->getDouble();
	}
	
	file << std::endl;
	file.close();
	
	file.open((outputFolder + obs + ".dat").c_str(), std::ios_base::out | std::ios_base::app);
	file << iteration << " " << readParams(obs)->getDouble() << std::endl; 
	file.close();
}
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
			throw std::runtime_error("Usage : CDMFT inputDirectory outputDirectory dataDirectory inputfilename iteration");
		
		std::string inputFolder = argv[1];
		std::string outputFolder = argv[2];
		std::string dataFolder = argv[3];
		std::string name = argv[4];
		int const iteration = std::atoi(argv[5]);
		std::string filename = "";
		std::string nodeName = "";
		if(iteration){
			filename = inputFolder + name + boost::lexical_cast<std::string>(iteration) + ".meas.json";
			nodeName = "Parameters";
		}else{
			filename = inputFolder + name + "0.json";
		}
		
		newIO::GenericReadFunc readParams(filename,nodeName);
		double const mu = readParams("mu")->getDouble();
		double const beta = readParams("beta")->getDouble();
		double const tpd = readParams("tpd")->getDouble();
		double const tpp = readParams("tpp")->getDouble();
		double tppp = tpp;
		//We want to read tppp if it is defined in the parameter file
		bool existstppp;
		const newIO::GenericReader* tpppRead = readParams("tppp",existstppp);
		if(tpppRead) {
			std::cout << "We have tppp different than tpp" << std::endl;
			tppp = tpppRead->getDouble();
		}
		//End of tppp read
		double const ep = readParams("ep")->getDouble();

		std::complex<double> w = .0;
		
		std::vector<RCuMatrix> selfEnergy;
		std::vector<RCuMatrix> hyb;
			

		//We read the Link file in order to know the structure of the self-energy and hybridation functions: 
		json_spirit::mArray jLink;
		readLinkFile(readParams("LINK")->getString(),jLink,outputFolder);
		//Now we create the map that will contain the components
		std::map<std::string,std::complex<double> > component_map;
		std::map<std::string,std::vector<std::pair<std::size_t,std::size_t> > > inverse_component_map;
		//First we create all different 
		std::size_t nSite_ = jLink.size()/2;
		for(std::size_t i=0;i<jLink.size();i++){
			for(std::size_t j=0;j<jLink.size();j++){
				component_map[jLink[i].get_array()[j].get_str()] = 0;
				if ( inverse_component_map.find(jLink[i].get_array()[j].get_str()) == inverse_component_map.end() ) {
					inverse_component_map[jLink[i].get_array()[j].get_str()] = std::vector<std::pair<std::size_t,std::size_t> >();
				}
				inverse_component_map[jLink[i].get_array()[j].get_str()].push_back(std::pair<std::size_t,std::size_t>(i,j));
			}
		}
		IO::WriteFunc writeHyb;

		if(iteration) {

			newIO::GenericReadFunc readMeas(inputFolder + name + boost::lexical_cast<std::string>(iteration) + ".meas.json","Measurements");
			readMeas.addSign(readMeas("Sign")->getDouble()); //Very important, otherwise the sign is not included in the simulation
			newIO::GenericReadFunc readHyb(outputFolder + boost::lexical_cast<std::string>(readParams("HYB")->getString()).c_str(),"");
			//We have to read all the Hyb components into variables so take them from LinkN.json
			//For that we need a table that stores the LinkN.json structure
			
			std::size_t NHyb = 0;
			std::size_t NGreen = 0;
			for (auto &p : component_map)
			{
				if(p.first != "empty"){ //We don't read the empty component
					if(NHyb == 0){
						NHyb = readHyb(p.first)->getSize();
						NGreen = readMeas("GreenI_" + p.first)->getSize();
					}else if(NHyb != readHyb(p.first)->getSize()){
						throw std::runtime_error(p.first + ": missmatch in entry length's of the hybridisation function.");
					}else if(NGreen != readMeas("GreenI_" + p.first)->getSize()){
						throw std::runtime_error(p.first + ": missmatch in entry length's of the measured Green's function function.");
					}
				}
			} 
			
			hyb.resize(NGreen);	
			//We initialize the Hybridization object from data		
			for(std::size_t n = 0; n < std::min(NHyb, NGreen); ++n) {
				//We iterate over all components and read them from the Hyb file
				for (auto &p : component_map)
				{
					if(p.first != "empty"){ //We don't read the empty component
						p.second = readHyb(p.first)->getFunction(n);
					}
				} 
				//We initialize the hybridization matrix according to the Link file.
				for(std::size_t i=0;i<jLink.size();i++){
					for(std::size_t j=0;j<jLink.size();j++){		
						std::complex<double> this_component = component_map[jLink[i].get_array()[j].get_str()];			
						hyb[n](i,j) = this_component;
						//We need to beware to initialize the down component according to the Nambu convention
						if(i >= nSite_ && j>= nSite_){
							hyb[n](i,j) = -std::conj(hyb[n](i,j));
						}
					}
				}				
			}
			//We initialize the Hyb object if the size doesn't match the Green's functions
			for (auto &p : component_map)
			{
				if(p.first != "empty"){ //We don't read the empty component
					p.second = readHyb(p.first)->getFM();
				}
			} 
			for(std::size_t n = std::min(NHyb, NGreen); n < NGreen; ++n) {
				
				std::complex<double> iomega(.0, M_PI*(2*n + 1)/beta);
				for(std::size_t i=0;i<jLink.size();i++){
					for(std::size_t j=0;j<jLink.size();j++){		
						std::complex<double> this_component = component_map[jLink[i].get_array()[j].get_str()]/iomega;		
						hyb[n](i,j) = this_component;
						//We need to beware to initialize the down component according to the Nambu convention
						if(i >= nSite_ && j>= nSite_){
							hyb[n](i,j) = -std::conj(hyb[n](i,j));
						}
					}
				}	
			}
			//We intialize the cluster Green's function
			std::vector<RCuMatrix> green(NGreen);
			for(std::size_t n = 0; n < NGreen; ++n) {
				for (auto &p : component_map)
				{
					if(p.first != "empty"){ //We don't read the empty component
						p.second = std::complex<double>(readMeas("GreenR_" + p.first)->getDouble(n),readMeas("GreenI_" + p.first)->getDouble(n));
					}
				} 
				for(std::size_t i=0;i<jLink.size();i++){
					for(std::size_t j=0;j<jLink.size();j++){		
						std::complex<double> this_component = component_map[jLink[i].get_array()[j].get_str()];			
						green[n](i,j) = this_component;
						//We need to beware to initialize the down component according to the Nambu convention
						if(i >= nSite_ && j>= nSite_){
							green[n](i,j) = -std::conj(green[n](i,j));
						}
					}
				}	
			}
			//We compute the selfEnergy
			for(std::size_t n = 0; n < NGreen; ++n) {
				std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
				
				RCuMatrix temp;
				temp(0, 0) = temp(1, 1) = temp(2, 2) = temp(3, 3) = iomega + mu;
				temp(4, 4) = temp(5, 5) = temp(6, 6) = temp(7, 7) = -std::conj(iomega + mu);

				temp -= hyb[n]; 	
				temp -= green[n].inv();
				
				selfEnergy.push_back(temp);
			}


			for (auto &p : inverse_component_map)
			{
				if(p.first != "empty"){
					writeHyb(p.first).FM() = readHyb(p.first)->getFM();
				}
			}


			//We read the observables into files
			bool existsn;
			const newIO::GenericReader* nRead = readMeas("n",existsn);
			if(nRead) {
				double const S = readMeas("S")->getDouble();
				readParams("mu")->setDouble(mu - S*(readMeas("N")->getDouble() - nRead->getDouble()));
			}
			
			{
				std::ofstream file(dataFolder + "sign.dat", std::ios_base::out | std::ios_base::app);
				file << iteration << " " << readMeas("Sign")->getDouble() << std::endl;
				file.close();
			}
			
			readScalSites("N", readMeas, iteration,dataFolder);
			readScalSites("k", readMeas, iteration,dataFolder);
			readScalSites("Sz", readMeas, iteration,dataFolder);
			readScalSites("D", readMeas, iteration,dataFolder);
			readScalSites("Chi0", readMeas, iteration,dataFolder);
						
			{
				
			    std::stringstream name; name << dataFolder << "pK" << iteration << ".dat";
			    std::ofstream file(name.str().c_str(), std::ios_base::out);
			    
				const newIO::GenericReader* pK_read = readMeas("pK");

			    for(unsigned int k = 0; k < pK_read->getSize(); ++k) 
			    	file << k << " " << pK_read->getDouble(k) << std::endl;
				
			    file.close();
			}
			
			if(readParams("EObs")->getInt() > .0) {
				
				{
					std::stringstream name; name << dataFolder << "ChiFullSites" << iteration << ".dat";
					std::ofstream file(name.str().c_str());
					
					for(unsigned int n = 0; n < readMeas("Chi")->getSize(); ++n) {
						file << 2*n*M_PI/beta;
						for(int i = 0; i < 4; ++i) {
							std::string s = boost::lexical_cast<std::string>(i);
							file << " " << readMeas("Chi_" + s)->getDouble(n);
						}
						file << std::endl;
					}
					
					file.close();
				}
				
				{
					std::stringstream name; name << dataFolder << "ChiFull" << iteration << ".dat";
					std::ofstream file(name.str().c_str());
				
					const newIO::GenericReader* Chi_read = readMeas("Chi");
				    for(unsigned int n = 0; n < Chi_read->getSize(); ++n){
				    		file << 2*n*M_PI/beta << " " << Chi_read->getDouble(n) << std::endl;
				    }
					file.close();
				}				
			};
			w = std::complex<double>(readParams("weightR")->getDouble(),readParams("weightI")->getDouble());
		} else {
			throw std::runtime_error("Please do not use this part of the program (iteration=0) as this file is not suited for starting a simulation.");
			//We initialize the simulation using a self file or an empty self-energy
			unsigned int const NSelf = beta*readParams("EGreen")->getInt()/(2*M_PI) + 1;
			
			std::ifstream selfFile(dataFolder + "self.dat");

			std::string dummy;
			selfEnergy.resize(NSelf);
			for(std::size_t n = 0; n < NSelf; ++n) {
				if(selfFile.good())  /* USING A COMPUTED ANORMAL SELFENERGY */
				{
					selfFile >> dummy;
					for (auto &p : component_map)
					{
						double real,imag;
						if(p.first != "empty"){ //We don't read the empty component
							selfFile >> real >> imag;
							p.second = std::complex<double>(real,imag);
						}
					} 
				}else /* USING A CUSTOM ANORMAL SELFENERGY */
				{
					double const delta = readParams("delta")->getDouble();
					double omega = (2*n + 1)*M_PI/beta;
					component_map["pphi"] = delta/(1. + omega*omega);
					component_map["mphi"] = -delta/(1. + omega*omega);
				}

				for(std::size_t i=0;i<jLink.size();i++){
					for(std::size_t j=0;j<jLink.size();j++){		
						std::complex<double> this_component = component_map[jLink[i].get_array()[j].get_str()];			

						selfEnergy[n](i,j) = this_component;
						//We need to beware to initialize the down component according to the Nambu convention
						if(i >= nSite_ && j>= nSite_){
							selfEnergy[n](i,j) = -std::conj(selfEnergy[n](i,j));
						}
					}
				}	
			}
			
			hyb.resize(selfEnergy.size());
			w = .0;			
		}

					
		std::ofstream selfFile((dataFolder + "self" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream greenFile((dataFolder + "green" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream hybFile((dataFolder + "hyb" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream structFile((dataFolder + "structure" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		
		Int::EulerMaclaurin2D<RCuMatrix> integrator(1.e-10, 4, 12);
		
		structFile << "Imaginary part of the Matsubara frequency";
		for (auto &p : inverse_component_map)
		{
			if(p.first != "empty"){
				structFile << " " << p.first;
			}
		}
		structFile << std::endl;
		for(std::size_t n = 0; n < selfEnergy.size(); ++n) {
			
			std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
			selfFile << iomega.imag();
			//We need to mean on all the indices included in the inverse component map
			//For each component type (00,01,11...) , we add the contributions of all the matrix coefficients that orrespond to those components
			for (auto &p : inverse_component_map)
			{
				
				std::complex<double> component_mean(0.,0.);
				std::size_t multiplicity = 0;
				if(p.first != "empty"){
					for(auto& pair : p.second){
						//Be careful of the Nambu convention
						if(pair.first >= nSite_ && pair.second >= nSite_){
							component_mean += -std::conj(selfEnergy[n](pair.first,pair.second));
						}else{
							component_mean += selfEnergy[n](pair.first,pair.second);
						}
						multiplicity+=1;
					}
					component_mean/=multiplicity;
					selfFile << " " << component_mean.real() << " " << component_mean.imag();
				}
			}
			selfFile << std::endl;
			
			RCuLatticeGreen latticeGreenRCu(iomega + mu, tpd, tpp, tppp, ep, selfEnergy[n]); 
			RCuMatrix greenNext = integrator(latticeGreenRCu, M_PI/2., M_PI/2.);

			greenFile << iomega.imag();
			//We do the same thing as for the selfEnergy above
			for (auto &p : inverse_component_map)
			{
				std::complex<double> component_mean(0.,0.);
				std::size_t multiplicity = 0;
				if(p.first != "empty"){
					for(auto& pair : p.second){
						//Be careful of the Nambu convention
						if(pair.first >= nSite_ && pair.second >= nSite_){
							component_mean += -std::conj(greenNext(pair.first,pair.second));
						}else{
							component_mean += greenNext(pair.first,pair.second);
						}
						multiplicity+=1;
					}
					component_mean/=multiplicity;
					greenFile << " " << component_mean.real() << " " << component_mean.imag();
				}
			}
			greenFile << std::endl;
			
			RCuMatrix hybNext;
			hybNext(0, 0) = hybNext(1, 1) = hybNext(2, 2) = hybNext(3, 3) = iomega + mu;
			hybNext(4, 4) = hybNext(5, 5) = hybNext(6, 6) = hybNext(7, 7) = -std::conj(iomega + mu);
			hybNext -= selfEnergy[n];
			hybNext -= greenNext.inv();

			hybFile << iomega.imag();
			//Finally we do that same procedure for the hyb file. It is a bit different this time as we also need to create th next Hyb"i".json file
			for (auto &p : inverse_component_map)
			{
				std::complex<double> component_mean_next(0.,0.);
				std::size_t multiplicity= 0;
				std::complex<double> component_mean_old(0.,0.);
				if(p.first != "empty"){
					for(auto& pair : p.second){
						//Be careful of the Nambu convention
						if(pair.first >= nSite_ && pair.second >= nSite_){
							component_mean_next += -std::conj(hybNext(pair.first,pair.second));
							component_mean_old += -std::conj(hyb[n](pair.first,pair.second));
						}else{
							component_mean_next += hybNext(pair.first,pair.second);
							component_mean_old += hyb[n](pair.first,pair.second);
						}
						multiplicity+=1;
					}
					component_mean_next/=multiplicity;
					component_mean_old/=multiplicity;
					hybFile << " " << component_mean_next.real() << " " << component_mean_next.imag();
					std::complex<double> toWriteNextIteration = (1. - w)*component_mean_next + w*component_mean_old;
					if(p.first == "pphi" || p.first == "mphi"){
						toWriteNextIteration = toWriteNextIteration.real();
					}
					writeHyb(p.first).push_back(toWriteNextIteration);
				}
			}
			hybFile << std::endl;
		}
		
		selfFile.close();
		greenFile.close();
		hybFile.close();
		//Old formula (kept for reference)
		/*
		writeHyb("00").FM() = 4*tpd*tpd; 
		writeHyb("01").FM() = -tpd*tpd; 
		writeHyb("11").FM() = .0; 
		*/


				
		readParams("HYB")->setString("Hyb" + boost::lexical_cast<std::string>(iteration + 1) + ".json");
		writeHyb.write(beta,outputFolder + "Hyb" + boost::lexical_cast<std::string>(iteration + 1) + ".json");

		readParams.write((outputFolder + "params" +  boost::lexical_cast<std::string>(iteration + 1) + ".json").c_str());
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
