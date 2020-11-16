#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"

/*Here are all the functions that will have to be replaced to get rid of ALPS 

readScal(name) (returns a Ut::MeasEntry, corresponding to the right string)
	For this function, we have to create a variable from json_spirit to store all those parameters, read from an xml file, we can replace that with a json file, no xml needed at all here (no alps)
read all parameters from the input params file (json file), alps is unneccesary here
readScal
readGreen
readVec

What would be interesting would be to keep approximatly the same structure of this file and to not use the IO.h anymore (no use)

*/
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
};


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
			std::cout << "We have tppp different than tpp" << std::endl;
			tppp = tpppRead->getDouble();
		}
		
		double const ep = readParams("ep")->getDouble();

		std::complex<double> w = .0;
		
		std::vector<RCuMatrix> selfEnergy;
		std::vector<RCuMatrix> hyb;
		if(iteration) {
			newIO::GenericReadFunc readMeas(inputFolder + name + boost::lexical_cast<std::string>(iteration) + ".meas.json","Measurements");
			readMeas.addSign(readMeas("Sign")->getDouble()); //Very important, otherwise the sign is not included in the simulation
			newIO::GenericReadFunc readHyb(outputFolder + boost::lexical_cast<std::string>(readParams("HYB")->getString()).c_str(),"");
			
			
			std::size_t const NHyb = readHyb("00")->getSize();
			if(NHyb != readHyb("01")->getSize()) throw std::runtime_error("01: missmatch in entry length's of the hybridisation function.");
			if(NHyb != readHyb("11")->getSize()) throw std::runtime_error("11: missmatch in entry length's of the hybridisation function.");
			
			std::size_t const NGreen = readMeas("GreenI_00")->getSize();
			
			hyb.resize(NGreen);			
			for(std::size_t n = 0; n < std::min(NHyb, NGreen); ++n) {
				std::complex<double> h00 = readHyb("00")->getFunction(n);
				std::complex<double> h01 = readHyb("01")->getFunction(n);
				std::complex<double> h11 = readHyb("11")->getFunction(n);
				
				hyb[n]("d_0", "d_0") = h00; hyb[n]("d_0", "d_1") = h01; hyb[n]("d_0", "d_2") = h11; hyb[n]("d_0", "d_3") = h01;
				hyb[n]("d_1", "d_0") = h01; hyb[n]("d_1", "d_1") = h00; hyb[n]("d_1", "d_2") = h01; hyb[n]("d_1", "d_3") = h11;
				hyb[n]("d_2", "d_0") = h11; hyb[n]("d_2", "d_1") = h01; hyb[n]("d_2", "d_2") = h00; hyb[n]("d_2", "d_3") = h01;
				hyb[n]("d_3", "d_0") = h01; hyb[n]("d_3", "d_1") = h11; hyb[n]("d_3", "d_2") = h01; hyb[n]("d_3", "d_3") = h00;
			}
			
			double FM00 = readHyb("00")->getFM();
			double FM01 = readHyb("01")->getFM();
			double FM11 = readHyb("11")->getFM();
			for(std::size_t n = std::min(NHyb, NGreen); n < NGreen; ++n) {
				
				std::complex<double> iomega(.0, M_PI*(2*n + 1)/beta);
				hyb[n]("d_0", "d_0") = FM00/iomega; hyb[n]("d_0", "d_1") = FM01/iomega; hyb[n]("d_0", "d_2") = FM11/iomega; hyb[n]("d_0", "d_3") = FM01/iomega;
				hyb[n]("d_1", "d_0") = FM01/iomega; hyb[n]("d_1", "d_1") = FM00/iomega; hyb[n]("d_1", "d_2") = FM01/iomega; hyb[n]("d_1", "d_3") = FM11/iomega;
				hyb[n]("d_2", "d_0") = FM11/iomega; hyb[n]("d_2", "d_1") = FM01/iomega; hyb[n]("d_2", "d_2") = FM00/iomega; hyb[n]("d_2", "d_3") = FM01/iomega;
				hyb[n]("d_3", "d_0") = FM01/iomega; hyb[n]("d_3", "d_1") = FM11/iomega; hyb[n]("d_3", "d_2") = FM01/iomega; hyb[n]("d_3", "d_3") = FM00/iomega;	
			}
			
			std::vector<RCuMatrix> green(NGreen);
			for(std::size_t n = 0; n < NGreen; ++n) {
				std::complex<double> g00 = std::complex<double>(readMeas("GreenR_00")->getDouble(n),readMeas("GreenI_00")->getDouble(n));
				std::complex<double> g01 = std::complex<double>(readMeas("GreenR_01")->getDouble(n),readMeas("GreenI_01")->getDouble(n));
				std::complex<double> g11 = std::complex<double>(readMeas("GreenR_11")->getDouble(n),readMeas("GreenI_11")->getDouble(n));
				
				green[n]("d_0", "d_0") = g00; green[n]("d_0", "d_1") = g01; green[n]("d_0", "d_2") = g11; green[n]("d_0", "d_3") = g01;
				green[n]("d_1", "d_0") = g01; green[n]("d_1", "d_1") = g00; green[n]("d_1", "d_2") = g01; green[n]("d_1", "d_3") = g11;
				green[n]("d_2", "d_0") = g11; green[n]("d_2", "d_1") = g01; green[n]("d_2", "d_2") = g00; green[n]("d_2", "d_3") = g01;
				green[n]("d_3", "d_0") = g01; green[n]("d_3", "d_1") = g11; green[n]("d_3", "d_2") = g01; green[n]("d_3", "d_3") = g00;				
			}

			IO::WriteFunc writeSelf;
			
			for(std::size_t n = 0; n < NGreen; ++n) {
				std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
				
				RCuMatrix temp = RCuMatrix::Diag(iomega + mu);
				temp -= hyb[n]; 	
				temp -= green[n].inv();
				
				selfEnergy.push_back(temp);
				
				writeSelf("00").push_back(temp("d_0", "d_0"));  
				writeSelf("01").push_back(temp("d_0", "d_1")); 
				writeSelf("11").push_back(temp("d_0", "d_2"));
			}
			
            writeSelf.write(beta, outputFolder + "Self" + boost::lexical_cast<std::string>(iteration) + ".json"); 
			

			bool existsn;
			const newIO::GenericReader* nRead = readMeas("n",existsn);
			if(nRead) {
				double const S = readMeas("S")->getDouble();
				const newIO::GenericReader* N_read = readMeas("N");
				readParams("mu")->setDouble(mu - S*(N_read->getDouble() - nRead->getDouble()));
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
			    for(unsigned int k = 0; k < pK_read->getSize(); ++k){
			    	file << k << " " << pK_read->getDouble(k) << std::endl;
			    }
			
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
			std::cout << w << std::endl;
		} else {
			selfEnergy.resize(beta*readParams("EGreen")->getInt()/(2*M_PI) + 1);
			hyb.resize(beta*readParams("EGreen")->getInt()/(2*M_PI) + 1);
			w = .0;
			
			/*
			IO::ReadFunc readSelf("Self.json");
			
			std::size_t const NSelf = readSelf("00").size();
			if(NSelf != readSelf("01").size()) throw std::runtime_error("01: missmatch in entry length's of the hybridisation function.");
			if(NSelf != readSelf("11").size()) throw std::runtime_error("11: missmatch in entry length's of the hybridisation function.");
			
			selfEnergy.resize(NSelf);
			hyb.resize(NSelf);
			
			for(std::size_t n = 0; n < NSelf; ++n) {
				std::complex<double> s00 = readSelf("00")[n];
				std::complex<double> s01 = readSelf("01")[n];
				std::complex<double> s11 = readSelf("11")[n];
				
				selfEnergy[n]("d_0", "d_0") = s00; selfEnergy[n]("d_0", "d_1") = s01; selfEnergy[n]("d_0", "d_2") = s11; selfEnergy[n]("d_0", "d_3") = s01;
				selfEnergy[n]("d_1", "d_0") = s01; selfEnergy[n]("d_1", "d_1") = s00; selfEnergy[n]("d_1", "d_2") = s01; selfEnergy[n]("d_1", "d_3") = s11;
				selfEnergy[n]("d_2", "d_0") = s11; selfEnergy[n]("d_2", "d_1") = s01; selfEnergy[n]("d_2", "d_2") = s00; selfEnergy[n]("d_2", "d_3") = s01;
				selfEnergy[n]("d_3", "d_0") = s01; selfEnergy[n]("d_3", "d_1") = s11; selfEnergy[n]("d_3", "d_2") = s01; selfEnergy[n]("d_3", "d_3") = s00;
			}
			
			w = .0;
			*/
		}
			//The real computations of the autocoherence relation stands right here
		{
			IO::WriteFunc writeHyb;
			
			std::ofstream selfFile((dataFolder + "self" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
			std::ofstream greenFile((dataFolder + "green" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
			std::ofstream hybFile((dataFolder + "hyb" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
			
			Int::EulerMaclaurin2D<RCuMatrix> integrator(1.e-10, 4, 12); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			for(std::size_t n = 0; n < selfEnergy.size(); ++n) {
				std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
				
				selfFile << iomega.imag() << " " 
				<< selfEnergy[n]("d_0", "d_0").real() << " " << selfEnergy[n]("d_0", "d_0").imag() << " " 
				<< selfEnergy[n]("d_0", "d_1").real() << " " << selfEnergy[n]("d_0", "d_1").imag() << " " 
				<< selfEnergy[n]("d_0", "d_2").real() << " " << selfEnergy[n]("d_0", "d_2").imag() << std::endl;
				
				/*Those two lines here compute the green function, the 2 most important lines of the program */
				//First we initialize the green's function with the right values
				RCuLatticeGreen latticeGreenRCu(iomega + mu, tpd, tpp, tppp, ep, selfEnergy[n]);
				//Then we compute the actual green's function, with the integral
				RCuMatrix greenNext = integrator(latticeGreenRCu, M_PI/2., M_PI/2.);
				
				greenFile << iomega.imag() << " "  
				<< greenNext("d_0", "d_0").real() << " " << greenNext("d_0", "d_0").imag() << " " 
				<< greenNext("d_0", "d_1").real() << " " << greenNext("d_0", "d_1").imag() << " "  
				<< greenNext("d_0", "d_2").real() << " " << greenNext("d_0", "d_2").imag() << std::endl;
				
				/*And you deduce te hybridization matrix from the green's function and the self energy */
				RCuMatrix hybNext = RCuMatrix::Diag(iomega + mu);
				hybNext -= selfEnergy[n];
				hybNext -= greenNext.inv();
				
				hybFile << iomega.imag() << " " 
				<< hybNext("d_0", "d_0").real() << " " << hybNext("d_0", "d_0").imag() << " " 
				<< hybNext("d_0", "d_1").real() << " " << hybNext("d_0", "d_1").imag() << " " 
				<< hybNext("d_0", "d_2").real() << " " << hybNext("d_0", "d_2").imag() << std::endl;
				
				writeHyb("00").push_back((1. - w)*hybNext("d_0", "d_0") + w*hyb[n]("d_0", "d_0")); 
				writeHyb("01").push_back((1. - w)*hybNext("d_0", "d_1") + w*hyb[n]("d_0", "d_1")); 
				writeHyb("11").push_back((1. - w)*hybNext("d_0", "d_2") + w*hyb[n]("d_0", "d_2"));
			}
			
			selfFile.close();
			greenFile.close();
			hybFile.close();
			
			writeHyb("00").FM() = 4*tpd*tpd; 
			writeHyb("01").FM() = -tpd*tpd; 
			writeHyb("11").FM() = .0; 

			readParams("HYB")->setString("Hyb" + boost::lexical_cast<std::string>(iteration + 1) + ".json");
			writeHyb.write(beta,outputFolder + "Hyb" + boost::lexical_cast<std::string>(iteration + 1) + ".json");
		}
			readParams.write((outputFolder + "params" +  boost::lexical_cast<std::string>(iteration + 1) + ".json").c_str());
	}
	catch(std::exception& exc) {
		std::cerr << exc.what() << "\n";
		std::cerr << "In CDMFT.cpp" << std::endl;
		return -1;
	}
	catch(...) {
		std::cerr << "Fatal Error: Unknown Exception!\n";
		return -2;
	}
	return 0;
}












