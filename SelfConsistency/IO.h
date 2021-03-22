#ifndef __NEWIO
#define __NEWIO

#include "nlohmann_json.hpp"
#include <fstream>
#include "Patrick/Utilities.h"
#include "Patrick/Plaquette/Plaquette.h"
using json=nlohmann::json;


bool exists(const json& j, const std::string& key)
{
    return j.find(key) != j.end();
}
/********/
/* Allows to divide all the elements in jObject by the value in jObject[key][0] */
/* jObject[key] keeps its value */
/* The jObject should be an object of arrays */
void divideAllBy(json& jObject, const std::string key){
	const double denominator = jObject[key][0];
	for (auto& el : jObject.items())
	{
		if(el.key() != key){
			/* Here el.value() is an array. We iterate over it to divide by the denominator variable*/
			for(size_t i=0;i<el.value().size();i++){
				double value = el.value()[i];
				value/=denominator;
				el.value()[i]=value;
			}
		}
	}
}
namespace IO{
	/*******************/
	/* LINK FILE UTILS */
	/*******************/
	/*********************/
	/* Loads a json file */
	/*********************/
	void readJsonFile(std::string fileName, json& jObject){
		std::cout << "Reading in " << fileName << std::endl;
	    std::ifstream file(fileName); 
	    if(file) {
	        file >> jObject;
	    }else{
	        throw std::runtime_error(fileName + " not found.");
	    }
	}
	/*********************/
	/*******************************/
	/* Write a json Object to file */
	/*******************************/
	void writeJsonToFile(std::string const fileName, json const& jObject){
		std::cout << "Writing to " << fileName << std::endl;
	    std::ofstream file(fileName); 
	    if(file) {
	        file << jObject.dump(1,'\t');
	    }else{
	        throw std::runtime_error("couldn't write to" + fileName);
	    }
	}
	/******************************/

	/***************************************/
	/****************************************/
	/* Write a Matsubara Observable to file */
	/****************************************/
	void writeMatsubaraToJsonFile(std::string const fileName, json& jObject,double const beta){
		for (auto& el : jObject.items()){
			el.value()["beta"] = beta;
			if(!exists(el.value(),"First Moment")){
				el.value()["First Moment"] = 0;
			}
			if(!exists(el.value(),"Second Moment")){
				el.value()["Second Moment"] = 0;
			}
		}
	    writeJsonToFile(fileName,jObject);
	}
	/******************************/
	/*************************************************************************************************************************/
	/* We read the Link file in order to know the structure of the self-energy and hybridation functions                     */
	/* In order to adapt to multiple simulation types, we check for different possible configurations                        */
	/*************************************************************************************************************************/
	void readLinkFromParams(json& jLink, std::string& outputFolder,json& jParams){
        /* First we check if there is the params file contains a LINK entry */
        /* If so, we read the whole Link file from there */
        if(exists(jParams,"LINK")){
        		std::string linkFileName = jParams["LINK"];
                readJsonFile(outputFolder + linkFileName,jLink);
        }else{

        /* Else, we assume the program provides a Normal and an Anomal part for the matrix structure. This form is spin symmetric */
        	std::string linkAFileName = jParams["LINKA"];
        	std::string linkNFileName = jParams["LINKN"];
            json jLinkA;
            json jLinkN;
            readJsonFile(outputFolder + linkAFileName,jLinkA);
            readJsonFile(outputFolder + linkNFileName,jLinkN);
            std::size_t nSite_ = jLinkN.size();
			jLink = jLinkN;            
            //First we create the full jLink matrix from the two little ones
            for(std::size_t i=0;i<nSite_;i++){
				jLink[i].insert(jLink[i].end(), jLinkA[i].begin(), jLinkA[i].end());				
				jLink.push_back(json::array());
				jLink[i+nSite_].insert(jLink[i+nSite_].end(), jLinkA[i].begin(), jLinkA[i].end());				
				jLink[i+nSite_].insert(jLink[i+nSite_].end(), jLinkN[i].begin(), jLinkN[i].end());				
            }
        }
	}
	/*************************************************************************************************************************/
	/****************************************************************************************/
	/* Distributes the components from component_map to matrix according to the jLink object */
	void component_map_to_matrix(json jLink,RCuMatrix& matrix,std::map<std::string,std::complex<double> >& component_map){
	    std::size_t nSite_ = jLink.size()/2;
	    for(std::size_t i=0;i<jLink.size();i++){
	        for(std::size_t j=0;j<jLink.size();j++){        
	            std::complex<double> this_component = component_map[jLink[i][j]];         
	            matrix(i,j) = this_component;
	            //We need to beware to initialize the down component according to the Nambu convention
	            if(i >= nSite_ && j>= nSite_){
	                matrix(i,j) = -std::conj(matrix(i,j));
	            }
	        }
	    }
	}
	/****************************************************************************************/
}
#endif
