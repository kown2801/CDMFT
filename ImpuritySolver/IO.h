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
    /* In order to adapt to multiple simulation types, we check for different possible configurations                        */
    /*************************************************************************************************************************/
    void readLinkFromParams(json& jLink, std::string& outputFolder,json const& jParams){
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
}
#endif