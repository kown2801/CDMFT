#ifndef __IO
#define __IO

#include "nlohmann_json.hpp"
#include <fstream>
#include "Patrick/Utilities.h"
#include "Patrick/Plaquette/Plaquette.h"
using json=nlohmann::json;

bool exists(const json& j, const std::string& key){
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
	/* JSON FILE UTILS */
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
	/*********************/
	/* Loads a json file at a specific iteration*/
	/*********************/
	/* From a Bundle file*/
	void readJsonFileInBundle(const std::string fileName, const int iteration, json & jObject){
		std::cout << "Reading in " << fileName << std::endl;
	    std::ifstream file(fileName); 
	    if(file) {
			std::string iterationString = std::to_string(iteration);
	        file >> jObject;
			jObject = jObject[iterationString];
			if(jObject.is_null()){
				throw std::runtime_error("Iteration " + std::to_string(iteration) + " not found in " + fileName);
			}
	    }else{
	        throw std::runtime_error(fileName + " not found.");
	    }
	}
	/* First from a plain file and then from a Bundle file*/
	void readJsonFileAtIteration(const std::string filename_start, const std::string filename_end, const int iteration, json & jObject){
		try{// We first check if the file exists on its own
			readJsonFile(filename_start + std::to_string(iteration) + filename_end,jObject);
		}catch( const std::runtime_error & e ){//Else we fetch it from a big `observableName` file
			std::cout << filename_start + std::to_string(iteration) + filename_end << " not found. Trying the bundle file." << std::endl;
			readJsonFileInBundle(filename_start + filename_end,iteration,jObject);
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
	
	/***************************************************************************/
	/* Saves the data in a big file (used to reduce the total number of files) */
	void writeInJsonDataFile(const std::string fileName, const int iteration, json const& jObject){
		//First we load the json file containing the already computed observableName from previous iterations (if there is any)
		json jExistingData;
		try{
			readJsonFile(fileName, jExistingData);
		}catch( const std::runtime_error & e ) {
                std::cout << fileName << " doesn't exist yet. We create it." << std::endl;
        }
		//Then we add the data from this iteration 
		std::string	iterationJsonIndex = std::to_string(iteration);
		jExistingData[iterationJsonIndex] = jObject;
		//And we save the whole object again in the file
		writeJsonToFile(fileName,jExistingData);
	}
	void writeInJsonDataFile(const std::string filename_start, const std::string filename_end, const int iteration, json const& jObject){
		writeInJsonDataFile(filename_start + filename_end,iteration,jObject);
	}
	/***************************************************************************/
	
	/****************************************/
	/* Write a Matsubara Observable to file */
	/****************************************/
	void formatMatsubaraData(json& jObject,double const beta){
		for (auto& el : jObject.items()){
			el.value()["beta"] = beta;
			if(!exists(el.value(),"First Moment")){
				el.value()["First Moment"] = 0;
			}
			if(!exists(el.value(),"Second Moment")){
				el.value()["Second Moment"] = 0;
			}
			el.value()["dataType"] = "FermionicMatsubaraFrequencies";
		}
	}
	/******************************/
	/*************************************************************************************************************************/
	/* We read the Link file in order to know the structure of the self-energy and hybridation functions                     */
	/*************************************************************************************************************************/
	void readLinkFromParams(json& jLink, std::string& outputFolder,json& jParams){
		if(exists(jParams,"LINK")){
			std::string linkFileName = jParams["LINK"];
			readJsonFile(outputFolder + linkFileName,jLink);
		}else{
			std::cout << "There is no LINK field in the parameter file" << std::endl;
			throw;
		}
	}
	/*************************************************************************************************************************/
	/****************************************************************************************/
	/* Distributes the components from component_map to matrix according to the jLink object */
	void component_map_to_matrix(json jLink,RCuMatrix& matrix,std::map<std::string,std::complex<double> >& component_map){
	    for(std::size_t i=0;i<jLink.size();i++){
	        for(std::size_t j=0;j<jLink.size();j++){        
	            std::complex<double> this_component = component_map[jLink[i][j]];         
	            matrix(i,j) = this_component;
	        }
	    }
	}
	/****************************************************************************************/
}
#endif
