#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"


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
		double const ep = readParams("ep")->getDouble();

		std::complex<double> w = .0;
		
		std::vector<RCuMatrix> selfEnergy;
		std::vector<RCuMatrix> hyb;
		if(iteration) {

			newIO::GenericReadFunc readMeas(inputFolder + name + boost::lexical_cast<std::string>(iteration) + ".meas.json","Measurements");
			readMeas.addSign(readMeas("Sign")->getDouble()); //Very important, otherwise the sign is not included in the simulation
			newIO::GenericReadFunc readHyb(outputFolder + boost::lexical_cast<std::string>(readParams("HYB")->getString()).c_str(),"");
			
			std::size_t const NHyb = readHyb("00")->getSize();
			if(NHyb != readHyb("01" )->getSize()) throw std::runtime_error("01: missmatch in entry length's of the hybridisation function.");
			if(NHyb != readHyb("11")->getSize()) throw std::runtime_error("11: missmatch in entry length's of the hybridisation function.");
			if(NHyb != readHyb("pphi")->getSize()) throw std::runtime_error("pphi: missmatch in entry length's of the hybridisation function.");
			if(NHyb != readHyb("mphi")->getSize()) throw std::runtime_error("mphi: missmatch in entry length's of the hybridisation function.");
			
			std::size_t const NGreen = readMeas("GreenI_00")->getSize();
			
			hyb.resize(NGreen);			
			for(std::size_t n = 0; n < std::min(NHyb, NGreen); ++n) {
				std::complex<double> h00 = readHyb("00")->getFunction(n);
				std::complex<double> h01 = readHyb("01")->getFunction(n);
				std::complex<double> h11 = readHyb("11")->getFunction(n);
				std::complex<double> pphi = readHyb("pphi")->getFunction(n);
				std::complex<double> mphi = readHyb("mphi")->getFunction(n);
				
				hyb[n]("d_0Up", "d_0Up") = h00; hyb[n]("d_0Up", "d_1Up") = h01; hyb[n]("d_0Up", "d_2Up") = h11; hyb[n]("d_0Up", "d_3Up") = h01;
				hyb[n]("d_1Up", "d_0Up") = h01; hyb[n]("d_1Up", "d_1Up") = h00; hyb[n]("d_1Up", "d_2Up") = h01; hyb[n]("d_1Up", "d_3Up") = h11;
				hyb[n]("d_2Up", "d_0Up") = h11; hyb[n]("d_2Up", "d_1Up") = h01; hyb[n]("d_2Up", "d_2Up") = h00; hyb[n]("d_2Up", "d_3Up") = h01;
				hyb[n]("d_3Up", "d_0Up") = h01; hyb[n]("d_3Up", "d_1Up") = h11; hyb[n]("d_3Up", "d_2Up") = h01; hyb[n]("d_3Up", "d_3Up") = h00;
				
				hyb[n]("d_0Down", "d_0Down") = -std::conj(h00); hyb[n]("d_0Down", "d_1Down") = -std::conj(h01); hyb[n]("d_0Down", "d_2Down") = -std::conj(h11); hyb[n]("d_0Down", "d_3Down") = -std::conj(h01);
				hyb[n]("d_1Down", "d_0Down") = -std::conj(h01); hyb[n]("d_1Down", "d_1Down") = -std::conj(h00); hyb[n]("d_1Down", "d_2Down") = -std::conj(h01); hyb[n]("d_1Down", "d_3Down") = -std::conj(h11);
				hyb[n]("d_2Down", "d_0Down") = -std::conj(h11); hyb[n]("d_2Down", "d_1Down") = -std::conj(h01); hyb[n]("d_2Down", "d_2Down") = -std::conj(h00); hyb[n]("d_2Down", "d_3Down") = -std::conj(h01);
				hyb[n]("d_3Down", "d_0Down") = -std::conj(h01); hyb[n]("d_3Down", "d_1Down") = -std::conj(h11); hyb[n]("d_3Down", "d_2Down") = -std::conj(h01); hyb[n]("d_3Down", "d_3Down") = -std::conj(h00);
				
				hyb[n]("d_0Up", "d_1Down") = 
				hyb[n]("d_0Down", "d_1Up") = 
				hyb[n]("d_1Up", "d_0Down") = 
				hyb[n]("d_1Down", "d_0Up") = 
				hyb[n]("d_2Up", "d_3Down") = 
				hyb[n]("d_2Down", "d_3Up") =
				hyb[n]("d_3Up", "d_2Down") = 
				hyb[n]("d_3Down", "d_2Up") = pphi;
				
				
				hyb[n]("d_1Up", "d_2Down") = 
				hyb[n]("d_1Down", "d_2Up") = 
				hyb[n]("d_2Up", "d_1Down") = 
				hyb[n]("d_2Down", "d_1Up") = 
				hyb[n]("d_3Up", "d_0Down") = 
				hyb[n]("d_3Down", "d_0Up") = 
				hyb[n]("d_0Up", "d_3Down") = 
				hyb[n]("d_0Down", "d_3Up") = mphi; 
				
			}
			
			for(std::size_t n = std::min(NHyb, NGreen); n < NGreen; ++n) {
				double FM00 = readHyb("00")->getFM();
				double FM01 = readHyb("01")->getFM();
				double FM11 = readHyb("11")->getFM();
				
				std::complex<double> iomega(.0, M_PI*(2*n + 1)/beta);
				hyb[n]("d_0Up", "d_0Up") = FM00/iomega; hyb[n]("d_0Up", "d_1Up") = FM01/iomega; hyb[n]("d_0Up", "d_2Up") = FM11/iomega; hyb[n]("d_0Up", "d_3Up") = FM01/iomega;
				hyb[n]("d_1Up", "d_0Up") = FM01/iomega; hyb[n]("d_1Up", "d_1Up") = FM00/iomega; hyb[n]("d_1Up", "d_2Up") = FM01/iomega; hyb[n]("d_1Up", "d_3Up") = FM11/iomega;
				hyb[n]("d_2Up", "d_0Up") = FM11/iomega; hyb[n]("d_2Up", "d_1Up") = FM01/iomega; hyb[n]("d_2Up", "d_2Up") = FM00/iomega; hyb[n]("d_2Up", "d_3Up") = FM01/iomega;
				hyb[n]("d_3Up", "d_0Up") = FM01/iomega; hyb[n]("d_3Up", "d_1Up") = FM11/iomega; hyb[n]("d_3Up", "d_2Up") = FM01/iomega; hyb[n]("d_3Up", "d_3Up") = FM00/iomega;
				
				hyb[n]("d_0Down", "d_0Down") = -std::conj(FM00/iomega); hyb[n]("d_0Down", "d_1Down") = -std::conj(FM01/iomega); hyb[n]("d_0Down", "d_2Down") = -std::conj(FM11/iomega); hyb[n]("d_0Down", "d_3Down") = -std::conj(FM01/iomega);
				hyb[n]("d_1Down", "d_0Down") = -std::conj(FM01/iomega); hyb[n]("d_1Down", "d_1Down") = -std::conj(FM00/iomega); hyb[n]("d_1Down", "d_2Down") = -std::conj(FM01/iomega); hyb[n]("d_1Down", "d_3Down") = -std::conj(FM11/iomega);
				hyb[n]("d_2Down", "d_0Down") = -std::conj(FM11/iomega); hyb[n]("d_2Down", "d_1Down") = -std::conj(FM01/iomega); hyb[n]("d_2Down", "d_2Down") = -std::conj(FM00/iomega); hyb[n]("d_2Down", "d_3Down") = -std::conj(FM01/iomega);
				hyb[n]("d_3Down", "d_0Down") = -std::conj(FM01/iomega); hyb[n]("d_3Down", "d_1Down") = -std::conj(FM11/iomega); hyb[n]("d_3Down", "d_2Down") = -std::conj(FM01/iomega); hyb[n]("d_3Down", "d_3Down") = -std::conj(FM00/iomega);
			}

			std::vector<RCuMatrix> green(NGreen);
			for(std::size_t n = 0; n < NGreen; ++n) {
				std::complex<double> g00 = std::complex<double>(readMeas("GreenR_00")->getDouble(n),readMeas("GreenI_00")->getDouble(n));
				std::complex<double> g01 = std::complex<double>(readMeas("GreenR_01")->getDouble(n),readMeas("GreenI_01")->getDouble(n));
				std::complex<double> g11 = std::complex<double>(readMeas("GreenR_11")->getDouble(n),readMeas("GreenI_11")->getDouble(n));
				std::complex<double> pphi = std::complex<double>(readMeas("GreenR_pphi")->getDouble(n),readMeas("GreenI_pphi")->getDouble(n));
				std::complex<double> mphi = std::complex<double>(readMeas("GreenR_mphi")->getDouble(n),readMeas("GreenI_mphi")->getDouble(n));
				
				
				green[n]("d_0Up", "d_0Up") = g00; green[n]("d_0Up", "d_1Up") = g01; green[n]("d_0Up", "d_2Up") = g11; green[n]("d_0Up", "d_3Up") = g01;
				green[n]("d_1Up", "d_0Up") = g01; green[n]("d_1Up", "d_1Up") = g00; green[n]("d_1Up", "d_2Up") = g01; green[n]("d_1Up", "d_3Up") = g11;
				green[n]("d_2Up", "d_0Up") = g11; green[n]("d_2Up", "d_1Up") = g01; green[n]("d_2Up", "d_2Up") = g00; green[n]("d_2Up", "d_3Up") = g01;
				green[n]("d_3Up", "d_0Up") = g01; green[n]("d_3Up", "d_1Up") = g11; green[n]("d_3Up", "d_2Up") = g01; green[n]("d_3Up", "d_3Up") = g00;	
				
				green[n]("d_0Down", "d_0Down") = -std::conj(g00); green[n]("d_0Down", "d_1Down") = -std::conj(g01); green[n]("d_0Down", "d_2Down") = -std::conj(g11); green[n]("d_0Down", "d_3Down") = -std::conj(g01);
				green[n]("d_1Down", "d_0Down") = -std::conj(g01); green[n]("d_1Down", "d_1Down") = -std::conj(g00); green[n]("d_1Down", "d_2Down") = -std::conj(g01); green[n]("d_1Down", "d_3Down") = -std::conj(g11);
				green[n]("d_2Down", "d_0Down") = -std::conj(g11); green[n]("d_2Down", "d_1Down") = -std::conj(g01); green[n]("d_2Down", "d_2Down") = -std::conj(g00); green[n]("d_2Down", "d_3Down") = -std::conj(g01);
				green[n]("d_3Down", "d_0Down") = -std::conj(g01); green[n]("d_3Down", "d_1Down") = -std::conj(g11); green[n]("d_3Down", "d_2Down") = -std::conj(g01); green[n]("d_3Down", "d_3Down") = -std::conj(g00);	
			
			    green[n]("d_0Up", "d_1Down") = 
				green[n]("d_0Down", "d_1Up") = 
				green[n]("d_1Up", "d_0Down") = 
				green[n]("d_1Down", "d_0Up") = 
				green[n]("d_2Up", "d_3Down") = 
				green[n]("d_2Down", "d_3Up") =
				green[n]("d_3Up", "d_2Down") = 
				green[n]("d_3Down", "d_2Up") = pphi;
				
				
				green[n]("d_1Up", "d_2Down") = 
				green[n]("d_1Down", "d_2Up") = 
				green[n]("d_2Up", "d_1Down") = 
				green[n]("d_2Down", "d_1Up") = 
				green[n]("d_3Up", "d_0Down") = 
				green[n]("d_3Down", "d_0Up") = 
				green[n]("d_0Up", "d_3Down") = 
				green[n]("d_0Down", "d_3Up") = mphi;
			}
			
			for(std::size_t n = 0; n < NGreen; ++n) {
				std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
				
				RCuMatrix temp;
				temp(0, 0) = temp(1, 1) = temp(2, 2) = temp(3, 3) = iomega + mu;
				temp(4, 4) = temp(5, 5) = temp(6, 6) = temp(7, 7) = -std::conj(iomega + mu);

				temp -= hyb[n]; 	
				temp -= green[n].inv();
				
				selfEnergy.push_back(temp);
			}

			
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
			/*
			unsigned int const NSelf = beta*static_cast<double>(params["EGreen"])/(2*M_PI) + 1;
			
			std::ifstream selfFile("self.dat");

			std::string dummy;
			selfEnergy.resize(NSelf);
			for(std::size_t n = 0; n < NSelf; ++n) {
				std::complex<double> s00;
				std::complex<double> s01;
				std::complex<double> s11;
				std::complex<double> pphi;
				std::complex<double> mphi;
				
				selfFile >> dummy >> s00.real() >> s00.imag() >> s01.real() >> s01.imag() >> s11.real() >> s11.imag() >> pphi.real() >> pphi.imag() >> mphi.real() >> mphi.imag();
				
				selfEnergy[n]("d_0Up", "d_0Up") = s00; selfEnergy[n]("d_0Up", "d_1Up") = s01; selfEnergy[n]("d_0Up", "d_2Up") = s11; selfEnergy[n]("d_0Up", "d_3Up") = s01;
				selfEnergy[n]("d_1Up", "d_0Up") = s01; selfEnergy[n]("d_1Up", "d_1Up") = s00; selfEnergy[n]("d_1Up", "d_2Up") = s01; selfEnergy[n]("d_1Up", "d_3Up") = s11;
				selfEnergy[n]("d_2Up", "d_0Up") = s11; selfEnergy[n]("d_2Up", "d_1Up") = s01; selfEnergy[n]("d_2Up", "d_2Up") = s00; selfEnergy[n]("d_2Up", "d_3Up") = s01;
				selfEnergy[n]("d_3Up", "d_0Up") = s01; selfEnergy[n]("d_3Up", "d_1Up") = s11; selfEnergy[n]("d_3Up", "d_2Up") = s01; selfEnergy[n]("d_3Up", "d_3Up") = s00;	
				
				selfEnergy[n]("d_0Down", "d_0Down") = -std::conj(s00); selfEnergy[n]("d_0Down", "d_1Down") = -std::conj(s01); selfEnergy[n]("d_0Down", "d_2Down") = -std::conj(s11); selfEnergy[n]("d_0Down", "d_3Down") = -std::conj(s01);
				selfEnergy[n]("d_1Down", "d_0Down") = -std::conj(s01); selfEnergy[n]("d_1Down", "d_1Down") = -std::conj(s00); selfEnergy[n]("d_1Down", "d_2Down") = -std::conj(s01); selfEnergy[n]("d_1Down", "d_3Down") = -std::conj(s11);
				selfEnergy[n]("d_2Down", "d_0Down") = -std::conj(s11); selfEnergy[n]("d_2Down", "d_1Down") = -std::conj(s01); selfEnergy[n]("d_2Down", "d_2Down") = -std::conj(s00); selfEnergy[n]("d_2Down", "d_3Down") = -std::conj(s01);
				selfEnergy[n]("d_3Down", "d_0Down") = -std::conj(s01); selfEnergy[n]("d_3Down", "d_1Down") = -std::conj(s11); selfEnergy[n]("d_3Down", "d_2Down") = -std::conj(s01); selfEnergy[n]("d_3Down", "d_3Down") = -std::conj(s00);	
				
				selfEnergy[n]("d_0Up", "d_1Down") = 
				selfEnergy[n]("d_0Down", "d_1Up") = 
				selfEnergy[n]("d_1Up", "d_0Down") = 
				selfEnergy[n]("d_1Down", "d_0Up") = 
				selfEnergy[n]("d_2Up", "d_3Down") = 
				selfEnergy[n]("d_2Down", "d_3Up") =
				selfEnergy[n]("d_3Up", "d_2Down") = 
				selfEnergy[n]("d_3Down", "d_2Up") = pphi;
				
				selfEnergy[n]("d_1Up", "d_2Down") = 
				selfEnergy[n]("d_1Down", "d_2Up") = 
				selfEnergy[n]("d_2Up", "d_1Down") = 
				selfEnergy[n]("d_2Down", "d_1Up") = 
				selfEnergy[n]("d_3Up", "d_0Down") = 
				selfEnergy[n]("d_3Down", "d_0Up") = 
				selfEnergy[n]("d_0Up", "d_3Down") = 
				selfEnergy[n]("d_0Down", "d_3Up") = mphi;
			}
			
			selfFile.close();
			/**/
			
			unsigned int const NSelf = beta*readParams("EGreen")->getInt()/(2*M_PI) + 1;
			
			std::ifstream selfFile(dataFolder + "self.dat");

			std::string dummy;
			selfEnergy.resize(NSelf);
			for(std::size_t n = 0; n < NSelf; ++n) {
				double s00_R, s00_I, s01_R, s01_I, s11_R, s11_I,pphi_R, pphi_I, mphi_R, mphi_I;
				if(selfFile.good())  /* USING A COMPUTED ANORMAL SELFENERGY */
				{
					selfFile >> dummy >> s00_R >> s00_I >> s01_R >> s01_I >> s11_R >> s11_I >> pphi_R >> pphi_I >> mphi_R >> mphi_I;
				}else /* USING A CUSTOM ANORMAL SELFENERGY */
				{
					double const delta = readParams("delta")->getDouble();
					double omega = (2*n + 1)*M_PI/beta;
					pphi_R = delta/(1. + omega*omega);
					mphi_R = -delta/(1. + omega*omega);
				}
					std::complex<double> s00(s00_R, s00_I);
					std::complex<double> s01(s01_R, s01_I);
					std::complex<double> s11(s11_R, s11_I);
					std::complex<double> pphi(pphi_R, pphi_I);
					std::complex<double> mphi(mphi_R, mphi_I);

			/**/
				selfEnergy[n]("d_0Up", "d_0Up") = s00; selfEnergy[n]("d_0Up", "d_1Up") = s01; selfEnergy[n]("d_0Up", "d_2Up") = s11; selfEnergy[n]("d_0Up", "d_3Up") = s01;
				selfEnergy[n]("d_1Up", "d_0Up") = s01; selfEnergy[n]("d_1Up", "d_1Up") = s00; selfEnergy[n]("d_1Up", "d_2Up") = s01; selfEnergy[n]("d_1Up", "d_3Up") = s11;
				selfEnergy[n]("d_2Up", "d_0Up") = s11; selfEnergy[n]("d_2Up", "d_1Up") = s01; selfEnergy[n]("d_2Up", "d_2Up") = s00; selfEnergy[n]("d_2Up", "d_3Up") = s01;
				selfEnergy[n]("d_3Up", "d_0Up") = s01; selfEnergy[n]("d_3Up", "d_1Up") = s11; selfEnergy[n]("d_3Up", "d_2Up") = s01; selfEnergy[n]("d_3Up", "d_3Up") = s00;	
				
				selfEnergy[n]("d_0Down", "d_0Down") = -std::conj(s00); selfEnergy[n]("d_0Down", "d_1Down") = -std::conj(s01); selfEnergy[n]("d_0Down", "d_2Down") = -std::conj(s11); selfEnergy[n]("d_0Down", "d_3Down") = -std::conj(s01);
				selfEnergy[n]("d_1Down", "d_0Down") = -std::conj(s01); selfEnergy[n]("d_1Down", "d_1Down") = -std::conj(s00); selfEnergy[n]("d_1Down", "d_2Down") = -std::conj(s01); selfEnergy[n]("d_1Down", "d_3Down") = -std::conj(s11);
				selfEnergy[n]("d_2Down", "d_0Down") = -std::conj(s11); selfEnergy[n]("d_2Down", "d_1Down") = -std::conj(s01); selfEnergy[n]("d_2Down", "d_2Down") = -std::conj(s00); selfEnergy[n]("d_2Down", "d_3Down") = -std::conj(s01);
				selfEnergy[n]("d_3Down", "d_0Down") = -std::conj(s01); selfEnergy[n]("d_3Down", "d_1Down") = -std::conj(s11); selfEnergy[n]("d_3Down", "d_2Down") = -std::conj(s01); selfEnergy[n]("d_3Down", "d_3Down") = -std::conj(s00);	

				selfEnergy[n]("d_0Up", "d_1Down") = 
				selfEnergy[n]("d_0Down", "d_1Up") = 
				selfEnergy[n]("d_1Up", "d_0Down") = 
				selfEnergy[n]("d_1Down", "d_0Up") = 
				selfEnergy[n]("d_2Up", "d_3Down") = 
				selfEnergy[n]("d_2Down", "d_3Up") =
				selfEnergy[n]("d_3Up", "d_2Down") = 
				selfEnergy[n]("d_3Down", "d_2Up") = pphi;
								
				selfEnergy[n]("d_1Up", "d_2Down") = 
				selfEnergy[n]("d_1Down", "d_2Up") = 
				selfEnergy[n]("d_2Up", "d_1Down") = 
				selfEnergy[n]("d_2Down", "d_1Up") = 
				selfEnergy[n]("d_3Up", "d_0Down") = 
				selfEnergy[n]("d_3Down", "d_0Up") = 
				selfEnergy[n]("d_0Up", "d_3Down") = 
				selfEnergy[n]("d_0Down", "d_3Up") = mphi;
			}
			
			hyb.resize(selfEnergy.size());
			w = .0;			
		}

		IO::WriteFunc writeHyb;
					
		std::ofstream selfFile((dataFolder + "self" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream greenFile((dataFolder + "green" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		std::ofstream hybFile((dataFolder + "hyb" + boost::lexical_cast<std::string>(iteration) + ".dat").c_str());
		
		Int::EulerMaclaurin2D<RCuMatrix> integrator(1.e-10, 4, 12);
		
		for(std::size_t n = 0; n < selfEnergy.size(); ++n) {
			
			std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);

			selfFile << iomega.imag() << " " 
			         << selfEnergy[n]("d_0Up", "d_0Up").real() << " " << selfEnergy[n]("d_0Up", "d_0Up").imag() << " " 
			         << selfEnergy[n]("d_0Up", "d_1Up").real() << " " << selfEnergy[n]("d_0Up", "d_1Up").imag() << " " 
			         << selfEnergy[n]("d_0Up", "d_2Up").real() << " " << selfEnergy[n]("d_0Up", "d_2Up").imag() << " "
			         << selfEnergy[n]("d_0Up", "d_1Down").real() << " " << selfEnergy[n]("d_0Up", "d_1Down").imag() << " "
			         << selfEnergy[n]("d_1Up", "d_2Down").real() << " " << selfEnergy[n]("d_1Up", "d_2Down").imag() << std::endl;
			
			RCuLatticeGreen latticeGreenRCu(iomega + mu, tpd, tpp, ep, selfEnergy[n]); 
			RCuMatrix greenNext = integrator(latticeGreenRCu, M_PI/2., M_PI/2.);

			greenFile << iomega.imag() << " "  
			          << greenNext("d_0Up", "d_0Up").real() << " " << greenNext("d_0Up", "d_0Up").imag() << " " 
			          << greenNext("d_0Up", "d_1Up").real() << " " << greenNext("d_0Up", "d_1Up").imag() << " " 
			          << greenNext("d_0Up", "d_2Up").real() << " " << greenNext("d_0Up", "d_2Up").imag() << " " 
			          << greenNext("d_0Up", "d_1Down").real() << " " << greenNext("d_0Up", "d_1Down").imag() << " "
			          << greenNext("d_1Up", "d_2Down").real() << " " << greenNext("d_1Up", "d_2Down").imag() << std::endl;
			
			RCuMatrix hybNext;
			hybNext(0, 0) = hybNext(1, 1) = hybNext(2, 2) = hybNext(3, 3) = iomega + mu;
			hybNext(4, 4) = hybNext(5, 5) = hybNext(6, 6) = hybNext(7, 7) = -std::conj(iomega + mu);
			hybNext -= selfEnergy[n];
			hybNext -= greenNext.inv();

			hybFile << iomega.imag() << " " 
			        << hybNext("d_0Up", "d_0Up").real() << " " << hybNext("d_0Up", "d_0Up").imag() << " " 
			        << hybNext("d_0Up", "d_1Up").real() << " " << hybNext("d_0Up", "d_1Up").imag() << " " 
			        << hybNext("d_0Up", "d_2Up").real() << " " << hybNext("d_0Up", "d_2Up").imag() << " "
			        << hybNext("d_0Up", "d_1Down").real() << " " << hybNext("d_0Up", "d_1Down").imag() << " "
			        << hybNext("d_1Up", "d_2Down").real() << " " << hybNext("d_1Up", "d_2Down").imag() << std::endl;
			
			writeHyb("00").push_back((1. - w)*hybNext("d_0Up", "d_0Up") + w*hyb[n]("d_0Up", "d_0Up")); 
			writeHyb("01").push_back((1. - w)*hybNext("d_0Up", "d_1Up") + w*hyb[n]("d_0Up", "d_1Up")); 
			writeHyb("11").push_back((1. - w)*hybNext("d_0Up", "d_2Up") + w*hyb[n]("d_0Up", "d_2Up"));
			writeHyb("pphi").push_back(((1. - w)*hybNext("d_0Up", "d_1Down") + w*hyb[n]("d_0Up", "d_1Down")).real());
			writeHyb("mphi").push_back(((1. - w)*hybNext("d_1Up", "d_2Down") + w*hyb[n]("d_1Up", "d_2Down")).real());
		}
		
		selfFile.close();
		greenFile.close();
		hybFile.close();

		writeHyb("00").FM() = 4*tpd*tpd; 
		writeHyb("01").FM() = -tpd*tpd; 
		writeHyb("11").FM() = .0; 
				
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












