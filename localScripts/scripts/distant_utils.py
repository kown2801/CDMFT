import select
import os
import numpy as np
from scripts import distant_consts as CONSTS

from io import BytesIO

scripts_folder = CONSTS.scripts_dir_from_main
all_data_folder = CONSTS.distant_data_dir
main_dir = CONSTS.distant_main_dir

#Reads a distant file into a file-like object (to prevent from creating temporary files into the directory)
def read_distant_file(c,file_path):
    io_obj = BytesIO()
    c.get(file_path,io_obj)
    io_obj.seek(0)
    return io_obj

#Loads the useful simulation data from the folder in the remote computer to here
#Returns a dict with some information about the simulation in folder
def get_folder_data(folder,c):
    folder_data = dict()
    folder_data["name"] = folder
    result = c.run(CONSTS.launch_numpy_actions + "-a order_parameter -f " + folder,hide=True)
    datafolder = os.path.join(all_data_folder,folder)
    folder_data["graph"] = np.load(read_distant_file(c,os.path.join(datafolder,"order_graph.npy")))
    folder_data["sign_graph"] = np.loadtxt(read_distant_file(c,os.path.join(datafolder,"DATA/sign.dat")))[:,1]
    folder_data["N_graph"] = np.loadtxt(read_distant_file(c,os.path.join(datafolder,"DATA/N.dat")))[:,1]
    try:
        folder_data["pn_graph"] = np.loadtxt(read_distant_file(c,os.path.join(datafolder,"DATA/pn.dat")))[:,1]
    except Exception as e:
        folder_data["pn_graph"] = None
        print(e)
        print("We couldn't read the pn file. Maybe it doesn't exist !")
    return folder_data

#Loading all the data files for each simulation 
#Returns nothing
#Adds to the all_folder_data array container one value per folder in all_folder_names
#If there was an error loading the data, the value is 0
#Otherwise the value is a dict with informations about the simulation
def get_all_folder_data(all_folder_data,all_folder_names,c):
    for f in all_folder_names:
        try:
            all_folder_data.append(get_folder_data(f,c))
        except Exception as e:
            print("Error at " + f + " : " + str(e))
            all_folder_data.append(0)
    print("Finished")
   