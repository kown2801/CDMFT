import numpy as np
import os
import json
import progressbar
import glob
import re
from AnalysisUtilities import utils
from AnalysisUtilities import constructG as c
from AnalysisUtilities import integrator


def find_missing_stiffness(input_dir,data_dir,convergedFor):
    try:
        lst = np.sort(np.loadtxt(os.path.join(data_dir,"stiffness.dat")).reshape(-1,2)[:,0].astype(int))
    except:
        lst = []
    should_be = utils.get_all_selfs(data_dir)
    max_iter = np.max(should_be)
    should_be = should_be[should_be > max_iter-convergedFor]
    to_compute = should_be[np.logical_not(np.isin(should_be,lst))]
    return to_compute

#This is the function that allows you to compute the stiffness for the convergedFor last iterations. This is done because it is long to compute the stiffness for all iterations locally. You should use a value here bigger (or equal) than the one used in Order-Parameter.ipynb in order to have sufficient computed iterations
def last_convergedFor_in_folder(files_dir,convergedFor):
    input_dir = os.path.join(files_dir,"IN/")
    output_dir = os.path.join(files_dir,"OUT/")
    data_dir = os.path.join(files_dir,"DATA/")
    all_data_dir = os.path.join(files_dir,"../")
    to_compute = find_missing_stiffness(input_dir,data_dir,convergedFor)
    if not to_compute.size == 0:
        print(files_dir)
        for el in progressbar.progressbar(to_compute):
            try:
                file = open(os.path.join(input_dir,"params" + str(el) + ".json"))
            except:
                file = open(os.path.join(input_dir,"params1.json"))
            params = json.load(file)
            file.close()
            stiffness = compute_stiffness(params,input_dir,os.path.join(data_dir,"self" + str(el) + ".json"))
            stiffness_filename = os.path.join(data_dir,"stiffness.dat")
            new_stiffness_data = np.array([el,stiffness]).reshape(-1,2)
            try:
                stiffness_data = np.loadtxt(stiffness_filename).reshape(-1,2)
                stiffness_data = np.append(stiffness_data,new_stiffness_data,axis=0)
                stiffness_data = stiffness_data[np.argsort(stiffness_data[:,0])]
            except Exception as e:
                stiffness_data = new_stiffness_data  
            np.savetxt(stiffness_filename,stiffness_data,fmt='%i %f')
    return 0


#Compute the stiffness for one self-energy file. This is done using the AnalysisUtilities/integrator.py file that is directly taken from Patrick Sémon's C++ integration code (EulerMaclaurin2D integration) and the AnalysisUtilities/constructG.py file that computes the terms of the integral in a parallel np.array fashion
def compute_stiffness(params,input_dir,file):
    if "tppp" not in params:
        params["tppp"] = 1
    integrand = c.Patrick_Integrand(params,input_dir,file)
    
    # Debut du decompte du temps
    error = 1e-2
    min_value = 1e-4
    stiffness = 0
    last_stiffness = 0
    z = 0
    while z < integrand.getLen() and (z == 0 or stiffness == 0 or last_stiffness/stiffness*(integrand.getLen() - z) > error/10) and (z==0 or stiffness>min_value) :
        integrand.setZ(z)
        #Because the sum is over all (positive and negative) Matsubara frequencies, there is a factor of 2 here.
        last_stiffness = 2*integrator.integrate(integrand.getStiffnessFunction(True),np.pi,np.pi,4,8,error)[0]
        stiffness += last_stiffness
        z+=1
    return stiffness