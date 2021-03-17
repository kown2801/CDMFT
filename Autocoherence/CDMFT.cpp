#include "Patrick/Integrators.h"
#include "newIO.h"
#include "Patrick/IO.h"
#include "Patrick/Plaquette/Plaquette.h"
#include <json_spirit.h>
    

/*****************************************************************************/
/* Alters the new Hyb file that comes out of the self-consistency relations. */
/* This is used to add constraints to the model (for example pphi and mphi are forced to be real)*/
std::complex<double> hyb_constraints(std::string component,std::complex<double> before){
    if(component == "pphi" || component == "mphi"){
        before = before.real();
    }
    return before;
}
/*****************************************************************************/
/**************************************************************************/
/* Initializes the selfEnergy when starting a new simulation from scratch */
/* This function may use all the parameter loaded in readParams           */
void initial_self_energy(newIO::GenericReadFunc& readParams,int n,std::map<std::string,std::complex<double> >& component_map){
    double const delta = readParams("delta")->getDouble();
    double omega = (2*n + 1)*M_PI/readParams("beta")->getDouble();
    component_map["pphi"] = delta/(1. + omega*omega);
    component_map["mphi"] = -delta/(1. + omega*omega);
}
/**************************************************************************/
/***************************************************************************/
/* Initializes the hyb moments when starting a new simulation from scratch */
/* This function may use all the parameter loaded in readParams            */
void initial_Hyb_moments(IO::WriteFunc& writeHyb,newIO::GenericReadFunc& readParams){
    double tpd = readParams("tpd")->getDouble();
    writeHyb("00").FM() = 4*tpd*tpd; 
    writeHyb("01").FM() = -tpd*tpd; 
    writeHyb("11").FM() = .0;
}
/***************************************************************************/
/****************************************************/
/* Post processing of the simple double observables */
void readScalSites(std::string obs, newIO::GenericReadFunc& readParams, int iteration,std::string outputFolder) {
    std::ofstream file((outputFolder + obs + "Sites.dat").c_str(), std::ios_base::out | std::ios_base::app);                
    file << iteration; 
    
    for(int i = 0; i < 4; ++i) {
        std::string s = std::to_string(i); 
        file << " " << readParams(obs + "_" + s)->getDouble();
    }
    
    file << std::endl;
    file.close();
    
    file.open((outputFolder + obs + ".dat").c_str(), std::ios_base::out | std::ios_base::app);
    file << iteration << " " << readParams(obs)->getDouble() << std::endl; 
    file.close();
}
/****************************************************/
/************************************************************************************/
/* Saves the data of matrix into the writeDat object that is used to save json data */
template<class T>
void write_Matsubara_data_to_file(IO::GenericWriteFunc<T>& writeDat, std::map<std::string,std::vector<std::pair<std::size_t,std::size_t> > >& inverse_component_map,RCuMatrix& matrix){
    //We need to mean on all the indices included in the inverse component map
    //For each component type (00,01,11...) , we add the contributions of all the matrix coefficients that orrespond to those components
    std::size_t nSite_ = matrix.size()/2;
    for (auto &p : inverse_component_map)
    {
        
        std::complex<double> component_mean(0.,0.);
        std::size_t multiplicity = 0;
        if(p.first != "empty"){
            for(auto& pair : p.second){
                //Be careful of the Nambu convention
                if(pair.first >= nSite_ && pair.second >= nSite_){
                    component_mean += -std::conj(matrix(pair.first,pair.second));
                }else{
                    component_mean += matrix(pair.first,pair.second);
                }
                multiplicity+=1;
            }
            component_mean/=multiplicity;
            writeDat(p.first).push_back(component_mean);
        }
    }
}
/************************************************************************************/

/***************************************************/
/* This scripts does the self-consistency relations*/
/* It also post-processes the observables, adding the sign and saving them files for future use */
int main(int argc, char** argv)
{
    try {
        if(argc != 6) 
            throw std::runtime_error("Usage : CDMFT inputDirectory outputDirectory dataDirectory inputfilename iteration");
        /******************************/
        /* Initialisation of variable */
        std::string inputFolder = argv[1];
        std::string outputFolder = argv[2];
        std::string dataFolder = argv[3];
        std::string name = argv[4];
        int const iteration = std::atoi(argv[5]);
        std::string filename = "";
        std::string nodeName = "";
        if(iteration){
            filename = inputFolder + name + std::to_string(iteration) + ".meas.json";
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
        /******************************/
            
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



        IO::WriteFunc writeThisHyb;
        IO::WriteFunc writeNextHyb;
        IO::WriteDat writeSelf;
        IO::WriteDat writeGreen;

        if(iteration) {

            newIO::GenericReadFunc readMeas(inputFolder + name + std::to_string(iteration) + ".meas.json","Measurements");
            readMeas.addSign(readMeas("Sign")->getDouble()); //Very important, otherwise the sign is not included in the simulation
            newIO::GenericReadFunc readHyb(outputFolder + readParams("HYB")->getString(),"");
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
            /************************************************/
            /**** We initialize the Hybridization object from data.  ******/
            for(std::size_t n = 0; n < std::min(NHyb, NGreen); ++n) {
                //We iterate over all components and read them from the Hyb file
                for (auto &p : component_map)
                {
                    if(p.first != "empty"){ //We don't read the empty component
                        p.second = readHyb(p.first)->getFunction(n);
                    }
                } 
                //We initialize the hybridization matrix according to the Link file.
                newIO::component_map_to_matrix(jLink,hyb[n],component_map);
            }
            /*** End initialisation from Matsubara data *****/
            /************************************************/
            /**************************************************************************************/
            /**** We initialize the Hyb object if the size doesn't match the Green's functions ****/
            /*** We use only the first moment expansion of the Hybridization for those components (this is usually a good starting point for the cycle) ***/
            for (auto &p : component_map)
            {
                if(p.first != "empty"){ //We don't read the empty component
                    p.second = readHyb(p.first)->getFM();
                }
            } 
            for(std::size_t n = std::min(NHyb, NGreen); n < NGreen; ++n) {
                std::complex<double> iomega(.0, M_PI*(2*n + 1)/beta);
                std::map<std::string,std::complex<double> > component_map_divided_by_iomega;
                for (auto &p : component_map)
                {
                    if(p.first != "empty"){ //We don't read the empty component
                        component_map_divided_by_iomega[p.first] = p.second/iomega;
                    }
                } 
                newIO::component_map_to_matrix(jLink,hyb[n],component_map_divided_by_iomega);
            }
            /*** End initialisation from Moments *****/
            /*****************************************/
            /********************************************************************************/
            /* We intialize the cluster Green's function from the Impurity Solver solution **/
            std::vector<RCuMatrix> green(NGreen);
            for(std::size_t n = 0; n < NGreen; ++n) {
                for (auto &p : component_map)
                {
                    if(p.first != "empty"){ //We don't read the empty component
                        p.second = std::complex<double>(readMeas("GreenR_" + p.first)->getDouble(n),readMeas("GreenI_" + p.first)->getDouble(n));
                    }
                } 
                newIO::component_map_to_matrix(jLink,green[n],component_map);
            }
            /* End Initialization of the cluster Green's function */
            /******************************************************/
            /*****************************/
            /* We compute the selfEnergy */
            for(std::size_t n = 0; n < NGreen; ++n) {
                std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
                
                RCuMatrix temp;
                for(std::size_t i = 0;i<nSite_;i++){
                    temp(i,i) = iomega + mu;
                    temp(i + nSite_,i + nSite_) = -std::conj(iomega + mu);
                }

                temp -= hyb[n];     
                temp -= green[n].inv();
                
                selfEnergy.push_back(temp);
            }
            /* End compute selfEnergy */
            /**************************/
            /**************************************/
            /* We read the observables into files */
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
                            std::string s = std::to_string(i);
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
            /* End Reading observables */
            /***************************/

            /*************************************************/
            /* We copy the first moment to the next Hyb file */
            for (auto &p : inverse_component_map)
            {
                if(p.first != "empty"){
                    writeNextHyb(p.first).FM() = readHyb(p.first)->getFM();
                }
            }
            /************************************************/


        } else {
            /******************************************************************************************************/
            /* Initialization of the self-energy using a self0.dat file or an the initialize_self_energy function */
            //We initialize the simulation using a self file or an empty self-energy
            unsigned int const NSelf = beta*readParams("EGreen")->getInt()/(2*M_PI) + 1;
            
            newIO::GenericReadFunc readSelf(dataFolder + "self0.json","");

            std::string dummy;
            selfEnergy.resize(NSelf);
            for(std::size_t n = 0; n < NSelf; ++n) {
                /*******************************************/
                /* We try loading the selfEnergy file */
                if(readSelf.good())  
                {
                    for (auto &p : component_map)
                    {
                        if(p.first != "empty"){ //We don't read the empty component
                            p.second = readSelf(p.first)->getFunction(n);
                        }
                    } 
                }else{
                    /*******************************************************************/
                    /* If this does not work, we initialize using this custom function */
                    initial_self_energy(readParams,n,component_map);
                    /*******************************************************************/
                }
                /*******************************************/
                newIO::component_map_to_matrix(jLink,selfEnergy[n],component_map);
            }
            
            hyb.resize(selfEnergy.size());
            initial_Hyb_moments(writeNextHyb,readParams);
            w = .0;         
        }
        /************************************************************************************************/
        /* Now we compute the next cluster Green's function and from that, the next hybridation function */
        Int::EulerMaclaurin2D<RCuMatrix> integrator(1.e-10, 4, 12);
        for(std::size_t n = 0; n < selfEnergy.size(); ++n) {
            
            std::complex<double> iomega(.0, (2*n + 1)*M_PI/beta);
        
            write_Matsubara_data_to_file(writeSelf,inverse_component_map,selfEnergy[n]);
            
            RCuLatticeGreen latticeGreenRCu(iomega + mu, tpd, tpp, tppp, ep, selfEnergy[n]); 
            RCuMatrix greenNext = integrator(latticeGreenRCu, M_PI/2., M_PI/2.);

            write_Matsubara_data_to_file(writeGreen,inverse_component_map,greenNext);
            
            RCuMatrix hybNext;
            for(std::size_t i = 0;i<nSite_;i++){
                hybNext(i,i) = iomega + mu;
                hybNext(i + nSite_,i + nSite_) = -std::conj(iomega + mu);
            }
            hybNext -= selfEnergy[n];
            hybNext -= greenNext.inv();
            /**********************************************************/
            /* Finally we save the data for the output Hyb.json file. */
            /* We need to do it once for the old Hyb file and one for the new in order to take a linear combination of both */
            write_Matsubara_data_to_file(writeThisHyb,inverse_component_map,hyb[n]);
            write_Matsubara_data_to_file(writeNextHyb,inverse_component_map,hybNext);
            /**********************************************************/
        }
        /************************************************************************************************/
        /************************************************************************/
        /*** Now we need to take the w weight into account for the next Hyb *****/
        w = std::complex<double>(readParams("weightR")->getDouble(),readParams("weightI")->getDouble());
        std::map<std::string, Hyb::Write>& hybNextMap = writeNextHyb.getMap();
        for (auto &p : hybNextMap){
            for(std::size_t i = 0;i<p.second.size();i++){
                std::complex<double> intermediate = (1. - w)*p.second[i] + w*writeThisHyb(p.first)[i];
                p.second[i] = hyb_constraints(p.first,intermediate);
            }
        }
        /************************************************************************/
        /**************************************/
        /* We save the self and green objects */
        writeSelf.write(beta,dataFolder + "self" + std::to_string(iteration) + ".json");
        writeGreen.write(beta,dataFolder + "green" + std::to_string(iteration) + ".json");
        /**************************************/

                
        /*******************************************************/
        /* We get ready for the next impurity Solver iteration */
        readParams("HYB")->setString("Hyb" + std::to_string(iteration + 1) + ".json");
        writeNextHyb.write(beta,outputFolder + "Hyb" + std::to_string(iteration + 1) + ".json");
        readParams.write(outputFolder + "params" +  std::to_string(iteration + 1) + ".json");
        /*******************************************************/
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









