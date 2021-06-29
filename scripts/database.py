import sqlite3
import json
from pathlib import Path
import numpy as np
from enum import Enum
from datetime import date
from importlib import reload
import bundle_utils as BU
import utils
import warnings

#Custom warnings to print only the message and not the lines or other informations
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return "Warning : " + str(msg) + '\n'

warnings.formatwarning = custom_formatwarning



#If you want to use this library, you should NEVER EVER move the folders around yourself. 
#The database works by using the relative paths of the simuation folders. Moving them breaks the database

single_valued_observables = ["Chi0","docc","ekin","k","n","np","sign","stiffness","Sz"]
site_valued_observables = ["Chi0Sites","doccSites","kSites","nSites","SzSites"]
observables_to_compute = ["occupation","superconducting_order_parameter"]
observables_not_in_database = ["green"]
sorting_parameters = ["dataset","hysteresis_direction"]

all_observables = single_valued_observables + site_valued_observables + observables_to_compute

default_values = {"tppp":1.0,"zero_order":0,"superconducting_order_parameter":0.0,"parameter_test":5.9875,"stiffness":0}



def compute_occupation(N_data,np_data):
    occupation = 2*N_data + 4*np_data
    return occupation

def compute_superconducting_order_parameter(G,parameters):
    beta = parameters["beta"]
    order = (G[:,0,:] - G[:,1,:]).real/2 * 2/beta
    order = np.sum(order,axis = 1)
    return order
    
#This dictionary variable *to_compute_data* is used to implement the on-the-fly calculations of quantities related to observables. 
#Those observables shouln't be matsubara frequency dependent
#To implement your own calculation, you need multiple things
    # First a name. Be careful not to choose twice the same name of something that exists already
    # A function that will be responsible for computing the quantity 
    # An array of the names of the observables you need in order to do your computation. 
    # In case you need the parameters for your computation indicate "Indices".
    # Be careful that an observable can't use the result of an other computed observable in its computation
    # For example here *superconducting_order_parameter* shouldn't be computed using *occupation* (even if it makes no senss)
    
    #The function takes for arguments the needed observables.
        # The remaining arguments will be array of observables that you specified in the "observables" field of the dict. 
        # For each observable, the first dimension is the iterations and the other dimensions are the dimension of the observables type (the copper occupation is just a number, the Green's function is an array along the matsubara frequencies and the different components)
    
#EXAMPLE :  the superconducting order parameter is computed from the anomal part of the Green's function 
#           and beta that is a parameter. So we indicate [{"green":["mphi","pphi"]},"Indices"]. 
#           The compute_order_parameter function has then three arguments. 
#           An array of all Green's function for all iterations and an dictionary of the simulation parameters.

to_compute_data = {"occupation" : {"compute" : compute_occupation, 
                                      "observables" :  ["n","np"]},
                   
                       "superconducting_order_parameter" : {"compute" : compute_superconducting_order_parameter,
                                  "observables" : [{"green":["mphi","pphi"]},"Indices"]
                                           }
                  }
     
#This allows to load the observable_name.dat files located in the dataFolder 
#It also corrects the order of the data in case the iterations are not in the right order 
def get_single_dat(dataFolder,observable_name,converged_from=None,max_iteration = None):
    data = np.loadtxt(dataFolder/(observable_name + ".dat"))
    indices = data[:,0]
    #We want to make sure the last converged_from iterations present in the file have indeed the right iteration indices
    if max_iteration != None:
        if converged_from == None:
            converged_from = max_iteration
        if not np.all(indices[-converged_from:] == np.arange(max_iteration-converged_from+1,max_iteration+1)):
            #try and sort the indices if possible and correct it in the file
            sorted_indices = np.argsort(indices)
            data = data[sorted_indices]
            indices = indices[sorted_indices]
            if not np.all(indices[-converged_from:] == np.arange(max_iteration-converged_from+1,max_iteration+1)):
                raise Exception(str(observable_name) + " misses some values in " + str(dataFolder) + ". The last " + str(converged_from) + " iterations do not match the self energy files\n","iteration_dont_match")
            else:
                np.savetxt(dataFolder/(observable_name + ".dat"),data,fmt="%i" + ' %f'*(data.shape[1]-1))
                warnings.warn("Corrected the order of iterations at " + str(observable_name) + " in " + str(dataFolder))
    obs_data = data[:,1:]  
    return obs_data    
    
#This is used to load complicated json file structure to arrays
def reduce_multiple_to_array(observable_data,observable_name,observable_structure,data_folder):
    #We iterate on the iterations
    is_numpy_array_data = False
    for i,el in enumerate(observable_data):
        array_data = []
        for j in observable_structure:
            if j in el:
                if "real" in el[j]:
                    is_numpy_array_data = True
                    array_data.append(np.array(el[j]["real"]) + 1j*np.array(el[j]["imag"]))
                else:
                    if isinstance(el[j],list):
                        array_data.append(np.array(el[j]))
                    else:
                        raise Exception("The type of data  at observable " + observable_name + " is not yet supported\
                        for extraction to a numpy array. We can't do the mean on it")                            
            else:
                raise Exception("No key " + str(j) + " for observable '" + str(observable_name)+ "' in folder " + str(data_folder))

            observable_data[i] = array_data
    if is_numpy_array_data:
        observable_data = np.array(observable_data)
    return observable_data

#This is a generic function that loads an observable depending on its datatype and storage type  
def load_observable(data_folder,observable_name,converged_from=None,max_iteration = None):
    if observable_name in single_valued_observables or observable_name in site_valued_observables :
        return get_single_dat(data_folder,observable_name,converged_from,max_iteration).squeeze()
    elif observable_name in observables_not_in_database:
        observable_data = BU.get_last_iterations_from_bundle(data_folder/observable_name,".json",converged_from,max_iteration)
        observable_structure = list(observable_data[0].keys())
        return reduce_multiple_to_array(observable_data,observable_name,observable_structure,data_folder)
    elif isinstance(observable_name,dict) and list(observable_name.keys())[0] in observables_not_in_database:
        observable_structure = observable_name
        observable_name = list(observable_name.keys())[0]
        observable_data = BU.get_last_iterations_from_bundle(data_folder/observable_name,".json",converged_from,max_iteration)
        return reduce_multiple_to_array(observable_data,observable_name,observable_structure[observable_name],data_folder)
    else:
        warnings.warn(observable_name  + " does not seem to exist")
        return None

# Computes the mean and the error on the last converged_from iterations for data (the first axis is always the iteration number)
def compute_mean_and_error(data,converged_from):
    data_mean = np.mean(data[-converged_from:],axis=0)
    data_error = np.std(np.abs(data[-converged_from:]),axis=0)/np.sqrt(converged_from)
    return data_mean,data_error
    
#This is a class to manage the simulation database  
class DatabaseManager:
    def __init__(self,database_file):
        self.connexion = sqlite3.connect(database_file)
        self.cursor = self.connexion.cursor()
        self.observable_fields = None
        self.folder_table_fields = None
        self.param_fields = None
        self.table_name = "simulation_folder"
        self.mean_table_name = "observable_mean"
        self.error_table_name = "observable_error"
        self.non_params_fields = ["date","folder"]
    
    def __del__(self):
        self.end()

    def end(self):
        try:
            # Save (commit) the changes
            self.commit()
            # We close the connection, because we are done with it
            self.connexion.close()
        except:
            pass
            
    #Save database changes
    def commit(self):
         self.connexion.commit()
         
    #Executes a query on the database        
    def execute(self,*args):
        return self.cursor.execute(*args)
        
    #Returns all the fields of the main (parameters) table
    def get_fields(self):
        if not self.folder_table_fields:
            this_cursor = self.execute("select * from " + self.table_name)
            self.folder_table_fields = list(map(lambda x: x[0], this_cursor.description))
        return self.folder_table_fields
    
    #Returns the parameters present in the database
    def get_params_fields(self):
        if not self.param_fields:
            self.param_fields = self.get_fields()
            for i in self.non_params_fields:
                try:
                    self.param_fields.remove(i)
                except:
                    pass
            for i in sorting_parameters:
                try:
                    self.param_fields.remove(i)
                except:
                    pass
        return self.param_fields
        
    #Returns the observables present in the database
    def get_observable_fields(self):
        if not self.observable_fields:
            this_cursor = self.execute("select * from " + self.mean_table_name)
            self.observable_fields = list(map(lambda x: x[0], this_cursor.description))
            self.observable_fields.remove("folder")
            self.observable_fields.remove("max_iter")
        return self.observable_fields
        
    #Verify that all registered parameters are in the database. Then, updates all data in the database
    def upgrade_parameters(self,param_type,converged_from=None):
        new_fields = set(param_type.keys()) - set(self.get_params_fields())
        #We add the mandatory sorting parameters
        for field in new_fields:
            if isinstance(param_type[field],int):
                datatype = "integer"
            elif isinstance(param_type[field],float):
                datatype = "real"
            else:
                datatype = "text"
            default_string = ""
            if field in default_values:
                default_string = "DEFAULT " + str(default_values[field]) + " NOT NULL"
            query = "ALTER TABLE " + self.table_name + " ADD COLUMN " + field + " " + datatype + " " + default_string
            self.execute(query)
        if converged_from != None:
            self.update_everything_in_table(converged_from)
        return new_fields
        
    #Verify that all registered observables are in the database. Then, updates all data in the database
    def upgrade_table(self,converged_from = None):
        new_fields = set(all_observables) - set(self.get_observable_fields())
        for field in new_fields:
            default_string = ""
            if field in default_values:
                default_string = "DEFAULT " + str(default_values[field]) + " NOT NULL"
            query = "ALTER TABLE " + self.mean_table_name + " ADD COLUMN " + field + " real " + default_string
            self.execute(query)
            query = "ALTER TABLE " + self.error_table_name + " ADD COLUMN " + field + " real " + default_string
            self.execute(query)
        if converged_from != None:
            self.update_everything_in_table(converged_from)
        return new_fields
        
    #Update all the registered folders
    def update_everything_in_table(self,converged_from):
        query = "SELECT folder FROM " + self.table_name
        query_result = list(self.execute(query))
        for folder in query_result:
            self.update(folder[0],converged_from)
            
    #Creates the 3 tables needed in the database : 
    #One for the parameters
    #One for the observable means
    #One for the observable errors
    def create_table(self,param_type):
        try:
            mysql_string = "CREATE TABLE " + self.table_name + " (date date, folder text, "
            #We add the mandatory sorting parameters
            mysql_string += " text, ".join(sorting_parameters)
            if len(sorting_parameters):
                mysql_string += " text"
            for i in param_type:
                if isinstance(param_type[i],int):
                    datatype = "integer"
                elif isinstance(param_type[i],float):
                    datatype = "real"
                else:
                    datatype = "text"
                default_string = ""
                if i in default_values:
                    default_string = "DEFAULT " + str(default_values[i]) + " NOT NULL"
                mysql_string += "," + i + " " + datatype + " " + default_string
            mysql_string += ",CONSTRAINT unique_folders PRIMARY KEY (folder))"
            self.execute(mysql_string)
        except sqlite3.OperationalError as e:
            warnings.warn("Table " + self.table_name + " exists already")

            
        common_mysql_string = " (folder text,max_iter integer,"\
                + ' real, '.join(single_valued_observables) + " real, "\
                + ' text, '.join(site_valued_observables) + " text, "\
                + ' real, '.join(observables_to_compute) + " real, "\
                + "CONSTRAINT unique_folders PRIMARY KEY (folder))"
        try:
            self.execute("CREATE TABLE " + self.mean_table_name + common_mysql_string)
        except sqlite3.OperationalError as e:
            warnings.warn("Table " + self.mean_table_name + " exists already")
        
        try:
            self.execute("CREATE TABLE " + self.error_table_name + common_mysql_string)
        except sqlite3.OperationalError as e:
            warnings.warn("Table " + self.error_table_name + " exists already")
    
    #Load the observables and parameters registered. This function is used to insert the observables in the database
    def load_and_format_observables(self,folder,converged_from,max_iter):
        data_folder = folder/"DATA"
        #First we prepare the parameters to store
        param_fields = self.get_params_fields()
        params = BU.get_in_bundle(folder/"IN"/"params",".json",1)
        params_to_insert = param_fields.copy()
        param_values = []
        for i in param_fields:
            if i in params:
                param_values.append(params[i])
            else:
                params_to_insert.remove(i)
                warnings.warn("When adding " + str(folder) + ", " + i + " was not found. Taking the default value")
        #Then we prepare the single-valued observables
        
        available_observables = all_observables.copy()
        observable_means = [] 
        observable_errors = [] 
        
        for observable_name in single_valued_observables:
            try:
                observable_data = load_observable(data_folder,observable_name,converged_from,max_iter)
                mean,error = compute_mean_and_error(observable_data,converged_from)
                observable_means.append(mean)
                observable_errors.append(error)
            except Exception as e:
                warnings.warn("Error for observable '" + observable_name + "' : " + str(e) + " , Taking the default value")
                available_observables.remove(observable_name)
                
            
        for observable_name in site_valued_observables:
            observable_data = load_observable(data_folder,observable_name,converged_from,max_iter)
            mean,error = compute_mean_and_error(observable_data,converged_from)
            observable_means.append(json.dumps(list(mean)))
            observable_errors.append(json.dumps(list(error)))
            
        for observable_name in observables_to_compute:
            if not observable_name in to_compute_data:
                warnings.warn("I don't know how to compute " + observable_name)
                return 
            #We try computing the observable
            try:
                #We gather the data needed for the computation

                args = []
                for i in to_compute_data[observable_name]["observables"]:
                    if i == "Indices":
                        args.append(params)
                    else:
                        args.append(load_observable(data_folder,i,converged_from,max_iter)) 


                #We compute the observable    
                observable_data = to_compute_data[observable_name]["compute"](*args)
                mean,error = compute_mean_and_error(observable_data,converged_from)
                observable_means.append(mean)
                observable_errors.append(error)
            except Exception as e:
                warnings.warn("Error : " + str(e))
                warnings.warn("Choosing the default value for observable " + observable_name + " as the data needed to compute it doesn't exist for folder " + str(folder))
                available_observables.remove(observable_name)
                
        return params_to_insert,param_values,available_observables,observable_means,observable_errors
        
    #Add a finished simulation folder to the database. This function is used to register a new folder
    def add_finished_simulation(self,folder,converged_from,sorting_parameters_values):
        folder = Path(folder).resolve()
        data_folder = folder/"DATA"
        max_iter = utils.get_max_iter(data_folder)
        #We make sure all the sorting_parameters_values are referenced
        sorting_parameters_query_part = ",".join(sorting_parameters) + ","*(len(sorting_parameters) != 0)
        sorting_parameters_query_values = []
        for i in sorting_parameters:
            if i not in sorting_parameters_values:
                raise Exception(i + " missing in the sorting_parameters_values argument while adding a simulation")
            else:
                sorting_parameters_query_values.append(sorting_parameters_values[i])
        
        params_to_insert,param_values,observable_to_insert,observable_means,observable_errors = self.load_and_format_observables(folder,converged_from,max_iter)
        try:
            a = self.execute("INSERT INTO " + self.table_name + " (date,folder," + sorting_parameters_query_part + \
                ",".join(params_to_insert) + ")" + " VALUES (" + ("?,")*(len(param_values)+len(sorting_parameters)+ 1) + "?)",
                          [date.today(),str(folder)] + sorting_parameters_query_values + param_values)
            self.execute("INSERT INTO " + self.mean_table_name + " (folder,max_iter," + ",".join(observable_to_insert) + ")"\
                      + " VALUES (" + ("?,")*(len(observable_means)+1) + "?)",[str(folder),str(max_iter)] + observable_means)
            self.execute("INSERT INTO " + self.error_table_name + " (folder,max_iter," + ",".join(observable_to_insert) + ")"\
                      + " VALUES (" + ("?,")*(len(observable_errors)+1) + "?)",[str(folder),str(max_iter)] + observable_errors)
        except sqlite3.IntegrityError as e:
            warnings.warn("Folder " + str(folder) + " already added to database")
     
    #Update an existing folder. It loads all registered parameters and observables from the simulation folder and updates accordingly 
    def update(self,folder,converged_from,sorting_parameters_values = {}):
        folder = Path(folder)
        #If the maximum iteration is the same or less than what we already have, we do nothing and print a message
        query_string = "SELECT max_iter FROM " + self.mean_table_name + " WHERE folder='" + str(folder) + "'"
        query_result = self.execute(query_string)
        try:
            old_max_iteration = query_result.fetchall()[0]
        except:
            warnings.warn("The folder " + str(folder) + " was not found in the database")
            return
        
        data_folder = folder/"DATA"
        current_max_iteration = utils.get_max_iter(data_folder)
        if old_max_iteration <= current_max_iteration:
            if old_max_iteration == current_max_iteration:
                warnings.warn("The folder " + str(folder) + " is already up-to-date. Updating anyway.")
            params_to_insert,param_values,observable_to_insert,observable_means,observable_errors = self.load_and_format_observables(folder,converged_from,current_max_iteration)
            #Now we check if we should update the sorting parameters for this simulation
            if sorting_parameters_values and (not set(sorting_parameters_values.keys()).issubset(set(sorting_parameters))):
                raise Exception("The sorting parameters were not recognized. " + str(set(sorting_parameters_values.keys()) - set(sorting_parameters)) + " doesn't exist in database")
            sorting_fields_in_db = []
            sorting_values_in_db = []
            for i in sorting_parameters_values:
                sorting_fields_in_db.append(i)
                sorting_values_in_db.append(sorting_parameters_values[i])
            try:
                query = "UPDATE " + self.table_name + " SET date = ?," + " = ?,".join(params_to_insert + sorting_fields_in_db) + " = ? WHERE folder='" + str(folder) + "'"
                self.execute(query,[date.today()] + param_values + sorting_values_in_db)
                query = "UPDATE " + self.mean_table_name + " SET " + " = ?,".join(observable_to_insert) + " = ? WHERE folder='" + str(folder) + "'"
                self.execute(query,observable_means)
                query = "UPDATE " + self.error_table_name + " SET " + " = ?,".join(observable_to_insert) + " = ? WHERE folder='" +str(folder) + "'"
                self.execute(query,observable_errors)
            except sqlite3.IntegrityError as e:
                warnings.warn("Update error on folder " + str(folder))
        else:
            warnings.warn("The folder " + str(folder) + " has less iterations than the entry in the database. Did you retry an point that was already computed ?")
    
    #Remove a folder from the database. It doesn't remove the files.
    def remove(self,folder):
        self.execute("DELETE FROM " + self.table_name + " WHERE folder='" + str(folder) + "'")
        self.execute("DELETE FROM " + self.mean_table_name + " WHERE folder='" + str(folder) + "'")
        self.execute("DELETE FROM " + self.error_table_name + " WHERE folder='" + str(folder) + "'")

    def get_params(self,conds):
        select_part = "SELECT * FROM " + self.table_name + " "
        query = select_part
        if conds:
            where_part = "WHERE "
            for i in conds:
                if isinstance(conds[i],str):
                    conds_text = "'" + conds[i] + "'"
                else:
                    conds_text  = str(conds[i])
                where_part += i + "=" + conds_text + " AND "
            where_part = where_part[:-5]
            query += where_part
        return self.execute(query)
    
    
    def get_that_there(self,to_plot,conds,converged_from = None):

        to_select = ""
        for i,el in enumerate(to_plot):
            if el in self.get_params_fields():
                to_select += "param." + el + ","
            elif el in self.get_observable_fields():
                to_select += "mean." + el + ",error." + el + ","
            elif el in observables_not_in_database or (isinstance(el,dict) and list(el.keys())[0] in observables_not_in_database):
                warnings.warn(str(el) + " not in database, but registered. We will load it directly from file.")
            else:
                raise ValueError(str(el) + " not available yet. Verify the spelling and the script file")
                
        to_select = to_select[:-1]
        select_part = "SELECT param.folder, mean.max_iter, " + to_select + " FROM " + self.table_name + " AS param"
        join_mean_part = "INNER JOIN " + self.mean_table_name + " AS mean ON mean.folder=param.folder"
        join_error_part = "INNER JOIN " + self.error_table_name + " AS error ON error.folder=param.folder"

        query = select_part + " " + join_mean_part + " " + join_error_part
        if conds:
            where_part = "WHERE "
            for i in conds:
                if isinstance(conds[i],str):
                    conds_text = "'" + conds[i] + "'"
                else:
                    conds_text  = str(conds[i])
                where_part += i + "=" + conds_text + " AND "
            where_part = where_part[:-5]
            query += " " + where_part
            
        query_result = self.execute(query)
        data_to_return = [np.empty((0,0))]*len(to_plot)
        for one_result in query_result:
            current_result_index = 2
            for index,el in enumerate(to_plot):
                if el in self.get_params_fields():
                    if data_to_return[index].shape[0] == 0:
                        data_to_return[index] = np.empty((0,2))
                    data_to_return[index] = np.vstack([data_to_return[index],[one_result[current_result_index],0]])
                    
                elif el in single_valued_observables or el in observables_to_compute:
                    if data_to_return[index].shape[0] == 0:
                        data_to_return[index] = np.empty((0,2))
                    data_to_return[index] = np.vstack([data_to_return[index],[one_result[current_result_index],one_result[current_result_index+1]]])
                    current_result_index+=1
                elif el in site_valued_observables:
                    mean_array = json.loads(one_result[current_result_index])
                    error_array = json.loads(one_result[current_result_index+1])
                    current_result_index += 1
                    
                    if data_to_return[index].shape[0] == 0:
                        data_to_return[index] = np.empty((0,len(mean_array),2))
                    this_cell = np.stack((mean_array,error_array)).transpose().reshape((1,-1,2))
                    data_to_return[index] = np.vstack([data_to_return[index],this_cell])
                
                current_result_index += 1
            for i,el in enumerate(to_plot):
                if el in observables_not_in_database or (isinstance(el,dict) and list(el.keys())[0] in observables_not_in_database):
                    if converged_from == None:
                        warnings.warn(el + " : in order to load data that doesn't fit in the db, you need to specify when the simulation has converged via the converged_from argument")
                        return
                    else:
                        observable_data = load_observable(Path(one_result[0])/"DATA",el,converged_from)
                        mean,error = compute_mean_and_error(observable_data,converged_from)
                        this_cell = np.stack((mean,error)).transpose().reshape((1,) + np.shape(mean) + (2,))
                        try:
                            if data_to_return[index].shape[0] == 0:
                                data_to_return[index] = np.empty((0,) + np.shape(this_cell)[1:])
                            data_to_return[index] = np.vstack([data_to_return[index],this_cell])
                        except ValueError as e:
                            warnings.warn(str(e))
                            warnings.warn("We couldn't load observable " + str(el) + " because the data has different shapes for different folders")
                
        return data_to_return