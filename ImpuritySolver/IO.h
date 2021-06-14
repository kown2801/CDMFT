#ifndef __NEWIO
#define __NEWIO

#include "nlohmann_json.hpp"
#include <fstream>
using json=nlohmann::json;


bool exists(const json& j, const std::string& key)
{
    return j.find(key) != j.end();
}

namespace IO
{
	/*********************/
    /* Loads a json file */
    /*********************/
    void readJsonFile(std::string fileName, json& jObject){
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
    void writeJsonToFile(std::string const fileName, json const & jObject,bool prettyPrint=true){
        std::ofstream file(fileName); 
        if(file) {
            if(prettyPrint){
                file << jObject.dump(1,'\t');
            }else{
                file << jObject.dump();
            }
        }else{
            throw std::runtime_error("couldn't write to" + fileName);
        }
    }
    /******************************/
	/*************************************************************************************************************************/
    /* We read the Link file in order to know the structure of the self-energy and hybridation functions                     */
    /*************************************************************************************************************************/
    void readLinkFromParams(json& jLink, std::string& outputFolder,json const& jParams){
        std::string linkFileName = jParams["LINK"];
		readJsonFile(outputFolder + linkFileName,jLink);
    }
    /*************************************************************************************************************************/
}
#endif