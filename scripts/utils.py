import itertools
import os
import json
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pickle
import bundle_utils as BU

data_save_file = "data-folders.back"

#Sorting functions helpers
def sort_points(E,N,key = None):
    E,N = zip(*sorted(zip(E, N),key=key))  
    return E,N

def sort_all(key=None,*args):
    return zip(*sorted(zip(*args),key=key)) 

#Returns all the data directory in the folder folder
#Uses a recursive function in order to search all folders inside *folder*
def find_all_parameters_in_folder(folder = "."):
    L = []
    def search_for_data(folder):
        current_folder = next(os.walk(folder))[0] 
        for i in next(os.walk(folder))[1]:
            try:
                folder_to_check = os.path.join(current_folder,i)
                if i in ["DATA","OUT","IN"]:
                    dataFolder = os.path.join(folder_to_check,"..","DATA")
                    json_data = BU.get_file_content(os.path.join(folder_to_check,"..","IN","params"),".json",1)
                    L.append([json_data,dataFolder])
                    return
                else:
                    search_for_data(folder_to_check)
            except Exception as e:
                print("Error at " + folder_to_check)
                print(str(e))
    search_for_data(folder)
    return L
 

#Retrieves all iteration numbers in dataFolder 
def get_all_selfs(data_folder):
    indices_str = BU.get_iterations_in_bundle(os.path.join(data_folder,"self"),".json")
    if indices_str == None: #Old version if no bundle is found in the folder
        all_selfs = glob.glob(os.path.join(data_folder,"self*.json"))
        all_indexes = []
        for i in all_selfs:
            all_indexes.append(int(re.findall(r'\d+', os.path.basename(i))[0]))
        return np.sort(np.array(all_indexes))
    return np.array(list(map(int, indices_str)))
    
#Retrieves the maximum iteration in dataFolder 
def get_max_iter(data_folder):
    all_selfs = get_all_selfs(data_folder)
    try:
        return(max(all_selfs))
    except Exception as e:
        print("There is no self files, take a look")
        print(str(e))

#Loads the components (list of strings) from a the measurement_name + ".json" file in dataFolder
#Rearranges it in a numpy array according to the order of the the components array. 
def load_matsubara_json(dataFolder,measurement_name,components,converged_from):
    max_iter = get_max_iter(dataFolder)
    i=1
    G = []
    if converged_from == None:
       converged_from = max_iter
    
    bundle_json_object = BU.get_bundle(os.path.join(dataFolder,measurement_name),".json")
    for i in range(max_iter-converged_from+1,max_iter+1):
        if bundle_json_object:
            json_object = bundle_json_object[str(i)]
        else:
            json_object = BU.get_file_content(os.path.join(dataFolder,measurement_name),".json",i)
        data = []
        for comp in components:
            if "real" in json_object[comp] :
                comp_values = np.array(json_object[comp]["real"]) + 1j*np.array(json_object[comp]["imag"])
            else:
                comp_values = np.array(json_object[comp])
            if len(data):
                comp_values = comp_values[..., np.newaxis]
                if len(data.shape) == 1:
                    data = data[..., np.newaxis]    
                data = np.concatenate((data,comp_values),axis=1)
            else:
                data = comp_values
        G.append(data)
    G = np.array(G,dtype=np.complex)
    return G
        
#Loads a matsubara .dat file.        
def load_matsubara_dat(dataFolder,measurement_name,converged_from):
    max_iter = get_max_iter(dataFolder)
    i=1
    G = []
    bundle_json_object = BU.get_bundle(os.path.join(dataFolder,measurement_name),".json")
    for i in range(max_iter-converged_from+1,max_iter+1):
        if bundle_json_object:
            json_object = bundle_json_object[str(i)]
        else:
            json_object = BU.get_file_content(os.path.join(dataFolder,measurement_name),".json",i)
        G_total = np.array(json_object)
        G.append(G_total)
    G = np.array(G)
    return G.squeeze()
