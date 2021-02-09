#ifndef __NEWIO
#define __NEWIO

#include <json_spirit.h>
#include <fstream>
namespace IO
{

	/*******************/
	/* LINK FILE UTILS */
	/*******************/
	/*********************/
	/* Loads a link file */
	void readLinkFile(std::string fileName, json_spirit::mArray& jLink,std::string inputFolder){
	    std::ifstream file(inputFolder + fileName); 
	    if(file) {
	        json_spirit::mValue temp;
	        json_spirit::read(file, temp); 
	        jLink = temp.get_array();
	    }else{
	        throw std::runtime_error(inputFolder + fileName + " not found.");
	    }
	}
	/*********************/
	/*************************************************************************************************************************/
	/* We read the Link file in order to know the structure of the self-energy and hybridation functions                     */
	/* In order to adapt to multiple simulation types, we check for different possible configurations                        */
	void readLinkFromParams(json_spirit::mArray& jLink, std::string& inputFolder, json_spirit::mObject const& jParams){
        /* First we check if there is the params file contains a LINK entry */
        /* If so, we read the whole Link file from there */
        if(jParams.find("LINK") != jParams.end()){
                readLinkFile(jParams.at("LINK").get_str(),jLink,inputFolder);
        }else{

        /* Else, we assume the program provides a Normal and an Anomal part for the matrix structure. This form is spin symmetric */
            json_spirit::mArray jLinkA;
            readLinkFile(jParams.at("LINKA").get_str(),jLinkA,inputFolder);
            json_spirit::mArray jLinkN;
            readLinkFile(jParams.at("LINKN").get_str(),jLinkN,inputFolder);
            
        	std::size_t nSite_ = jLinkN.size();
            jLink = json_spirit::mArray(2*nSite_);
            //First we create the full jLink matrix from the two little ones
            for(std::size_t i=0;i<nSite_;i++){
            jLink[i] = json_spirit::mArray(2*nSite_);
                jLink[i + nSite_] = json_spirit::mArray(2*nSite_);
                for(std::size_t j=0;j<nSite_;j++){
                    jLink[i].get_array()[j] = jLinkN[i].get_array()[j].get_str();
                    jLink[i + nSite_].get_array()[j + nSite_] = jLinkN[i].get_array()[j].get_str();
                    jLink[i + nSite_].get_array()[j] = jLinkA[i].get_array()[j].get_str();
                    jLink[i].get_array()[j + nSite_] = jLinkA[i].get_array()[j].get_str();
                }
            }
        }
	}
}
/*************************************************************************************************************************/
#endif