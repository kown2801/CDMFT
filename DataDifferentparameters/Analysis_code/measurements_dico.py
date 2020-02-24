import collections
import utils 

def get_all_matching(indexes,conds):
    matching = np.array([True for i in indexes])
    for key,value in conds.items():
        matching &= np.squeeze(indexes[key] == value)
    return matching

def to_list(object1):
    X,Y = object1
    return list(X) 
 
def get_order_parameter(green,beta):
    #We have data of the form green(iwn)(site)
    return (sum((green[:,4] - green[:,5])/2)).real * 2 / beta 

def get_plottable_order_parameter(greens,betas):
    order_parameters = []
    for i,item in enumerate(greens):
        order_parameters.append(get_order_parameter(item,betas[i]))       
    return order_parameters
        
class Measurements_Dico:
    def __init__(self,measurement_list):
        self.dict = {}
        for i in measurement_list["Normal"]:
            self.dict[i] = []
        for i in measurement_list["Matsubara"]:
            self.dict[i] = []
        self.dict["Indexes"] = []
        self.dict["Measurements_name"] = measurement_list
        self.plottable = 0
        self.sens_plottable = []
        self.param_type = 0
        
    def add(self,data,parameters = ""):
        #add the data to the dictionary 
        for i,measurement in data.items():
            self.dict[i].append(measurement)

        #create the data structure for parameters 
        if not self.param_type:
            self.param_type = ()
            for key,value in parameters.items():
                self.param_type += ((key,"float32"),)
            self.param_type = list(self.param_type)
            self.dict["Indexes"] = np.zeros(0,dtype=self.param_type)
            
        x = np.zeros(1,dtype=self.param_type)
        for key,value in parameters.items():
            try:
                x[key] = value
            except:pass
        self.dict["Indexes"] = np.append(self.dict["Indexes"],x)

    def is_computed_function(self,function):
        return function in self.dict["Measurements_name"]["Normal"] or function in self.dict["Measurements_name"]["Matsubara"]
    
    def plot_dico(self, function, axeX, orderBy = None, conds = {}):
        after_function = 0
        accepted_functions = ["order"]
        associated_functions = ["green"]
        
        if not self.is_computed_function(function):
            #We asked for a function of measurements
            if(function in accepted_functions and self.is_computed_function(associated_functions[accepted_functions.index(function)])):
                after_function = function
                function = associated_functions[accepted_functions.index(after_function)]
            else:
                raise ValueError("The function asked does not exist. Only the measurement in the list or 'order' are accepted")
        good_indexes = get_all_matching(self.dict["Indexes"],conds)
        selected_data = np.array(self.dict[function])[good_indexes]
        all_indexes = self.dict["Indexes"][good_indexes]
        if not orderBy:
            indexes = all_indexes[axeX]
            assert len(np.unique(indexes)) == len(indexes), "Duplicate values in indexes. \n" + str(indexes)
            indexes, selected_data = np.array(utils.sort_points(indexes,selected_data))
        else:
            indexes = all_indexes[[orderBy,axeX]]
            #The to_list function is used to get an ordered plot for the data
            indexes,selected_data = np.array(utils.sort_points(indexes,selected_data,to_list))
        
        if after_function :#We want to apply a function on the result
            if after_function == "order":
                if not orderBy:
                    e,betas = utils.sort_points(all_indexes[axeX],all_indexes["beta"])
                else:
                    e,betas = utils.sort_points(all_indexes[[orderBy,axeX]],all_indexes["beta"],to_list)
                selected_data = get_plottable_order_parameter(selected_data,betas)
       
        return indexes,selected_data
    
    
def get_plottable_arrays(X,Y):
    A,B = zip(*X)
    a = np.array([A,B]).transpose()
    b = np.cumsum(np.unique(a[:, 0], return_counts=True)[1])[:-1]
    indexes = np.split(a, b)
    selected_data = np.split(Y, b)
    return indexes,selected_data
    
 
import numpy as np
import utils
import greenperiodized as gp

def get_data_from_dict(dataFolder,measurement_list,converged_from):
    measurement_results = {}
    for observable_name in measurement_list["Normal"]:
        measurement_results[observable_name] = np.mean(np.loadtxt(dataFolder + observable_name + ".dat")[:,1][converged_from:])
    for observable_name in measurement_list["Matsubara"]:
        measurement_results[observable_name] = gp.load_matsubara_measurements_from_simulations(dataFolder,observable_name,converged_from)
    return measurement_results

def get_measurements_along(measurement_list,directory,converged_from):
    measurements_dico = Measurements_Dico(measurement_list)
    L = utils.find_all_parameters_in_folder(directory)
    for i in L:
        data = get_data_from_dict(i[1],measurement_list,converged_from)
        measurements_dico.add(data,i[0]["Parameters"])
    return measurements_dico