import select
import os
import numpy as np
from scripts import distant_consts as CONSTS
scripts_folder = CONSTS.scripts_dir_from_main
all_data_folder = CONSTS.distant_data_dir
main_dir = CONSTS.distant_main_dir

def get_folder_data(folder,c):
    folder_data = dict()
    folder_data["nom"] = folder
    result = c.run(CONSTS.launch_numpy_actions + "-a order_parameter -f " + folder,hide=True)
    datafolder = os.path.join(all_data_folder,folder)
    remote_file = c.get(os.path.join(datafolder,"order_graph.npy"))
    folder_data["graph"] = np.load("order_graph.npy")
    remote_sign_file = c.get(os.path.join(datafolder,"DATA/sign.dat"))
    folder_data["sign_graph"] = np.loadtxt("sign.dat")[:,1]
    remote_N_file = c.get(os.path.join(datafolder,"DATA/N.dat"))
    folder_data["N_graph"] = np.loadtxt("N.dat")[:,1]
    try:
        remote_pn_file = c.get(os.path.join(datafolder,"DATA/pn.dat"))
        folder_data["pn_graph"] = np.loadtxt("pn.dat")[:,1]
    except Exception as e:
        print(e)
        print("We couldn't read the pn file. Maybe it doesn't exist !")
    return folder_data

#Loading all the data files for each simulation (folder)
def get_all_folder_data(all_folder_data,all_folder_names,c):
    for f in all_folder_names:
        try:
            all_folder_data.append(get_folder_data(f,c))
        except Exception as e:
            print("Error at " + f + " : " + str(e))
            all_folder_data.append(0)
    print("Finished")
   