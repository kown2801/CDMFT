import database
from pathlib import Path
import bundle_utils as BU
import utils
import os
import shutil
import pickle
import json
import warnings

database_folder = Path.home()/"long_term_storage"/"AllData"
database_file = database_folder/"database.db"

folder_organisation = ["dataset","hysteresis_direction","ep","U","tpd","tppp","beta"]
#Sorting parameter example : 
#sorting_parameters_values = {"dataset":"only_supra","hysteresis_direction":"increasing"}

def find_available_data_dir_name(basename):
    current_name = str(basename) + ""
    i = 1
    while Path(current_name).is_dir():
        current_name = str(basename) + "_" + str(i)
        i+=1
    return Path(current_name)   
    
def move_and_tidy_folder(folder,all_data_folder,sorting_parameters_values):
    folder = Path(folder)
    all_data_folder = Path(all_data_folder)
    if "dataset" not in sorting_parameters_values and "hysteresis_direction" not in sorting_parameters_values:
        raise Exception("You need to indicate 'dataset' and 'hystersis_direction' in sorting_parameters_values if you want to move the folder")
    try:
        params = BU.get_in_bundle(folder/"IN"/"params",".json",1)
        folder_path = all_data_folder
        for el in folder_organisation:
            if el in sorting_parameters_values:
                folder_path/=str(sorting_parameters_values[el])
            else:
                if el in params:
                    folder_path/=(str(el) + "_" + str(params[el]))
                else:
                    folder_path/=(str(el) + "_" + str(database.default_values[el]))
    except Exception as e:
        print("Error: " + str(e) + " ... Proceeding with the default directory")
        folder_path = all_data_folder/"Default_folder"
    
    folder_path.mkdir(exist_ok=True,parents=True)
    if not folder == (folder_path/folder.name):
        folder_path/=find_available_data_dir_name(folder_path/folder.name).name
        #We need shutil for a cross platform move
        shutil.move(folder,folder_path)
    else:
         folder_path = folder
    return folder_path




def create_database():
    #Create the database using a parameter file
    param_type = BU.get_in_bundle(database_folder/"all_params",".json",0)
    c = database.DatabaseManager(database_file)
    c.create_table(param_type)
    c.end()
   
def add_simulation_folder(folder,converged_from,sorting_parameters_values):
    #Add all simulations folders from a folder in the database
    c = database.DatabaseManager(database_file)
    c.add_finished_simulation(folder,converged_from,sorting_parameters_values)
    c.end()
    
def add_all_in_folder(folder,converged_from,sorting_parameters_values):
    #Add all simulations folders from a folder in the database
    all_folders = utils.find_all_parameters_in_folder(folder)
    c = database.DatabaseManager(database_file)
    for i in all_folders:
        c.add_finished_simulation(Path(os.path.normpath(Path(i[1])/"..")),converged_from,sorting_parameters_values)
        c.commit()
    c.end()
    
def move_and_add_folder(folder,converged_from,sorting_parameters_values):
    c = database.DatabaseManager(database_file)
    new_folder_name = move_and_tidy_folder(folder,database_folder,sorting_parameters_values)
    c.add_finished_simulation(new_folder_name,converged_from,sorting_parameters_values)
    c.end() 
    
def move_and_add_all_in_folder(folder,converged_from,sorting_parameters_values):
    all_folders = utils.find_all_parameters_in_folder(folder)
    c = database.DatabaseManager(database_file)
    for i in all_folders:
        new_folder_name = move_and_tidy_folder(Path(os.path.normpath(Path(i[1])/"..")),database_folder,sorting_parameters_values)
        c.add_finished_simulation(new_folder_name,converged_from,sorting_parameters_values)
        c.commit()
    c.end() 
    
    
def get_all_simulations(conds = {}):
    #Get everything (in raw form) from the folder table in the db with conditions
    c = database.DatabaseManager(database_file)
    results = list(c.get_params(conds))
    print("There are " + str(len(results)) + " entries in the database" + " for those conditions"*(conds==False))
    c.end()
    return results

def get_connection():
    return database.DatabaseManager(database_file)
    
def get_nb_simulations(conds={}):
    c = database.DatabaseManager(database_file)
    results = list(c.get_params(conds))
    c.end()
    return len(results)
    
def remove(folder):
    c = database.DatabaseManager(database_file)
    c.remove(folder)
    c.end()
    
#Two Helper functions to sort the data
def takeFirstError(tuples):
    return [tuples[i][0] for i in range(0,len(tuples)-1)]
def sort_all(key=None,*args):
    return zip(*sorted(zip(*args),key=key))

#Gets the data in database, sorts it and sends it back to the local computer
def get_that_there_to_local(pickled_arguments):
    to_plot,conds,converged_from,orderBy = pickle.loads(str.encode(pickled_arguments))
    to_plot_local = to_plot.copy()
    if orderBy != None and orderBy not in to_plot:
        to_plot_local.append(orderBy)
    c = database.DatabaseManager(database_file)
    result = c.get_that_there(to_plot_local,conds,converged_from)
    c.end()
    for i in range(len(result)):
        result[i] = result[i].tolist()
    #Now we sort the data according to the orderBy parameter
    if orderBy != None and orderBy not in to_plot:
        orderBy_data = result[-1]
        to_plot_data = result[:-1]  
    else:
        if orderBy == None:
            orderBy = to_plot_local[0]
        orderBy_data = result[to_plot_local.index(orderBy)]
        to_plot_data = result
    to_plot_data = list(sort_all(takeFirstError,orderBy_data,*to_plot_data))[1:]
    #And we send back the data
    warnings.warn("Sending the data back")
    print(json.dumps(pickle.dumps(to_plot_data, protocol=0).decode('unicode_escape'))) # protocol 0 is printable ASCII
