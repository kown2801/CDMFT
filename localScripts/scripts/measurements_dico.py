import collections
from scripts import utils 
from scripts import greenperiodized as gp
import sys
import numpy as np
import os
import json
from copy import deepcopy

def compute_occupation(dataFolder,N_data,pn_data):
    occupation = 2*N_data + 4*pn_data
    save_to_dat(occupation,os.path.join(dataFolder,"occupation.dat"))

def compute_order_parameter(dataFolder,G,parameters):
    beta = parameters["beta"]
    order = (G[:,:,3] - G[:,:,4]).real/2 * 2/beta
    order = - np.sum(order,axis = 1)
    save_to_dat(order,os.path.join(dataFolder,"order.dat"))

accepted_functions = {"occupation" : {"compute" : compute_occupation, 
                                      "observables" :  {"Matsubara":[],"Normal" : ["N","pn"]}
                                     },
                       "order" : {"compute" : compute_order_parameter,
                                  "observables" : {"Matsubara":["green"],"Normal" : ["Indices"]}
                                }
                     }
#This dictionary variable *accepted_functions* is used to implement the on-the-fly calculations of quantities related to observables. Those observables shouln't be matsubara frequency dependent
#To implement your own calculation, you need multiple things
    #First a name. Be careful not to choose twice the same name of something that exists already
    # A function that takes for arguments the dataFolder and other observables.
        # dataFolder is a string that contains the DATA folder of the current simulation. The remaining arguments will be array of observables. The first dimension is the iterations and then the other dimensions are the dimension of the observables type (the copper occupation is just a number, the Green's function is an array along the matsubara frequencies and the different components). This function should compute the observable for each iteration and then save it to a observable_name.dat file in the dataFolder in the usual observable format "iteration data1 data2 ..." (see for example N.dat or Chi0Sites.dat )
    # An array of the names of the observables you need in order to do your computation. In case you need the parameters for your computation indicate "Indices". Be careful that you load all the right observables for your computation when creating the Measurements_dico or you will get an error. The exception is for Matsubara observables. As those are large arrays, if you don't need them anywhere else, you don't need to load them. This will be done for you. 

#EXAMPLE : the order parameter is computed from the Green's function and beta that is a parameter. So we indicate ["green","Indices"]. The compute_order_parameter function has then three two arguments. The dataFolder string, an array of all Green's function for all iterations and an array of the simulation parameters.



#Returns an array of boolean of the same shape as *indexes* that indicates if the elements correspond to *conds*
def get_all_matching(indexes,conds):
    matching = np.array([True for i in indexes])
    for key,value in conds.items():
        matching &= np.squeeze(np.logical_or(\
            np.logical_and(key == "tppp" and value==1.,indexes[key] == 0),\
            indexes[key] == value))
    return matching

def takeFirstError(tuples):
    return [tuples[i][0] for i in range(0,len(tuples)-1)]

def takeFirst(tuples):
    return [tuples[i] for i in range(0,len(tuples)-1)]

def save_to_dat(data,filename):
    n = len(data)
    data_to_file = np.zeros((n,2))
    data_to_file[:,0] = np.linspace(1,n,n).astype(int)
    data_to_file[:,1] = data
    np.savetxt(filename,data_to_file,fmt='%i %f')

class Measurements_Dico:
    def __init__(self,measurement_list):
        self.errors = False
        self.dict = {}
        for i in measurement_list["Normal"]:
            self.dict[i] = []
        for i in measurement_list["Matsubara"]:
            self.dict[i] = []
        self.dict["Indices"] = []
        self.dict["Measurements_name"] = measurement_list
        self.plottable = 0
        self.sens_plottable = []
        self.param_type = 0

#This function adds a point to the dictionary. It treats all the asked measurements but also all the parameters of the simulation        
    def add(self,data,parameters = ""):
        #add the data to the dictionary 
        for i,measurement in data.items():
            try:
                self.dict[i].append(measurement)
            except:
                pass

        #create the data structure for parameters 
        if not self.param_type:
            self.param_type = ()
            if("tppp" not in parameters):
                parameters["tppp"] = 1
            for key,value in parameters.items():
                self.param_type += ((key,"float64"),)
            self.param_type = list(self.param_type)
            self.dict["Indices"] = np.zeros(0,dtype=self.param_type)
            
        x = np.zeros(1,dtype=self.param_type)
        for key,value in parameters.items():
            try:
                x[key] = value
            except:pass
        self.dict["Indices"] = np.append(self.dict["Indices"],x)

    def is_computed_function(self,function):
        return function in self.dict["Measurements_name"]["Normal"] or function in self.dict["Measurements_name"]["Matsubara"]
    
    def get_function(self,good_indices,function):
        selected_indices = self.dict["Indices"][good_indices]
        if function in selected_indices.dtype.names:#The function asked for is a parameter
            to_return = selected_indices[function]
            if self.errors:
                to_return = np.expand_dims(to_return,axis=1)
                to_return = np.append(to_return,np.zeros((len(to_return),1)),axis=1)
            return to_return 
        elif self.is_computed_function(function):#The function asked for is a measurement
            return np.array(self.dict[function])[good_indices]
        else: #The function asked for is a non_computed quantity
            raise Exception(str(function) + " is not computed yet. Please include it in the measurements variable when loading the dictionary. If you want to implement something that does not already exist, please see in scripts/measurements_dico.py.")

#This function retrieves the data and puts it in the right form for it to be plotted in 2D or 3D with or without errors
#Ordering the points by the axeX or by another quantity using orderBy
#You can restrict the plot using the conds (ep=9 for example is a condition)
#You can retrieve the errors also
                                                                                  
    def plot_dico(self, to_plot,orderBy = None, conds = {}):
        #We always order by the first axis if not specified
        if orderBy == None:
            orderBy = to_plot[0]
                
        #Local copy of to_plot because lists are passed by reference
        to_plot_local = to_plot.copy()
        to_plot_local.append(orderBy)
        
        #We get the simulations corresponding to conds
        good_indexes = get_all_matching(self.dict["Indices"],conds)
        
        #We get the data corresponding to those simulations
        to_plot_data = []
        for axis in to_plot_local:
            to_plot_data.append(self.get_function(good_indexes,axis))
    
        #We sort according to orderBy data (the sorting function changes if the error is computed or not)
        orderBy_data = to_plot_data[-1]
        to_plot_data = to_plot_data[:-1]
        if self.errors:
            to_plot_data = list(utils.sort_all(takeFirstError,orderBy_data,*to_plot_data))[1:]
        else:
            to_plot_data = list(utils.sort_all(takeFirst,orderBy_data,*to_plot_data))[1:]
            
        #We put the data in the right shape for use
        for i,j in enumerate(to_plot_data):
            to_plot_data[i] = np.array(j).squeeze()                 
        return to_plot_data

#This allows to load the observable_name.dat files located in the dataFolder 
#and take the average over the last converged_from points
def get_single_dat(dataFolder,observable_name):
    obs_data = np.loadtxt(os.path.join(dataFolder,observable_name + ".dat"))[:,1:]  
    return obs_data

#This function is used to retrieve the quantities from a single directory "dataFolder/../". 
#parameters is the dict containing the parameters for the simulation
#measurements_list is the dict containing all the measurements we want to retrieve
#The simulation is supposed to be converged on the last converged_from iterations
#errors tells if we want to compute the error from the data.
def get_data_from_dict(dataFolder,parameters,measurement_list,converged_from,errors=False):
    measurement_results = {}
    compute_when_all_loaded = [] #Variable to compute quantities when all other quantities are loaded
    measurement_list_local = deepcopy(measurement_list)
    for observable_name in measurement_list_local["Normal"]:
        try:
            measurement_results[observable_name] = get_single_dat(dataFolder,observable_name).squeeze()
        except Exception as e:#There is an error reading the observable, maybe it is because we have to compute it from other 
            if observable_name in accepted_functions:
                for i in accepted_functions[observable_name]["observables"]["Matsubara"]:
                    measurement_list_local["Matsubara"].append(i)
                for i in accepted_functions[observable_name]["observables"]["Normal"]:
                    if i not in measurement_list_local["Normal"] and i != "Indices":
                        raise Exception("You need to include " + i + " in the normal observables in order to compute " + observable_name)
                compute_when_all_loaded.append(observable_name)
            elif observable_name=="stiffness":
                print("No stiffness yet for " + dataFolder)
                if errors:
                    measurement_results[observable_name] = [-1,0]*converged_from
                else:
                    measurement_results[observable_name] = [-1]*converged_from
    for observable_name in measurement_list_local["Matsubara"]:
        try:
            G = gp.load_matsubara_measurements_from_simulations(dataFolder,observable_name,converged_from)
        except Exception as e: #There is an error reading a Matsubara file so we report it
            print(str(e)  + " in " + dataFolder + " for observable " + observable_name)
        if observable_name in measurement_list:
            G_mean = np.mean(G[-converged_from:],axis=0)
            if errors: #We adapt in order to include the error on the array
                H = np.std(G[-converged_from:].real,axis=0)/np.sqrt(converged_from)
                H = np.expand_dims(H,axis=-1)
                G_mean = np.expand_dims(G_mean,axis=-1)
                G_mean = np.append(G_mean,H,axis=-1)
        else:
            G_mean = G
        measurement_results[observable_name] = G_mean
    for observable_name in compute_when_all_loaded:# Now we can compute observables that are computed from other
        args = []
        for i in accepted_functions[observable_name]["observables"]["Matsubara"]:
            args.append(measurement_results[i])
        for i in accepted_functions[observable_name]["observables"]["Normal"]:
            if i == "Indices":
                args.append(parameters)
            else:
                args.append(measurement_results[i])
        try:
            accepted_functions[observable_name]["compute"](dataFolder,*args)
        except Exception as e:
            raise Exception(str(e) + " in trying to compute " + observable_name + " from " + str(accepted_functions[observable_name]["observables"]) + ". I guess that the computed data don't have the right shapes.")
        measurement_results[observable_name] = get_single_dat(dataFolder,observable_name).squeeze()
    #After loading all the useful measurements, we compute the mean and the error if needed
    for k, v in measurement_results.items():
        v_mean = np.mean(v[-converged_from:],axis=0)
        if errors: #We adapt in order to include the error on the array
            H = np.std(np.abs(v[-converged_from:]),axis=0)/np.sqrt(converged_from)
            try:
                H = np.expand_dims(H,axis=-1)
                v = np.expand_dims(v_mean,axis=-1)
                v = np.append(v,H,axis=-1)
            except:
                v = np.array([v,H])
        else:
             v = v_mean
        measurement_results[k] = v
    return measurement_results


#This function is used to load the data you need to plot. It searches the directory to find all the data and loads
#all the functions in the mesaurement_list variable. 
#As we are dealing with DMFT results, converged_from corresponds to the number of iterations we keep to average on 
def get_measurements_along(measurement_list,directory,converged_from,errors=False):
    measurements_dico = Measurements_Dico(measurement_list)
    L = utils.find_all_parameters_in_folder(directory)
    for i in L:
        try:
            data = get_data_from_dict(i[1],i[0],measurement_list,converged_from,errors)
        except Exception as e:
            print("Error on folder " + i[1] + " : ")
            print(str(e))
        measurements_dico.add(data,i[0])
    measurements_dico.errors = errors
    return measurements_dico