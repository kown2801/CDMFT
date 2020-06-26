import itertools
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pickle
data_save_file = "data-folders.back"

#Returns all the data directory in the folder folder
#Uses a recursive function in order to search all folders inside *folder*
def find_all_parameters_in_folder(folder = "."):
    L = []
    def search_for_data(folder):
        current_folder = next(os.walk(folder))[0] 
        for i in next(os.walk(folder))[1]:
            try:
                folder_to_check = current_folder + "/" + i
                if i in ["DATA","OUT","IN"]:
                    dataFolder = os.path.join(folder_to_check,"../DATA/")
                    f = open(os.path.join(folder_to_check,"../IN/params1.json"),"r")
                    json_data = json.loads(f.read())
                    L.append([json_data,dataFolder])
                    return
                else:
                    search_for_data(folder_to_check)
            except:
                print("Error at " + folder_to_check)
    search_for_data(folder)
    return L
        


def sort_points(E,N,key = None):
    E,N = zip(*sorted(zip(E, N),key=key))  
    return E,N

def sort_all(key=None,*args):
    return zip(*sorted(zip(*args),key=key)) 

def contourPlot_BZ(kx_array, ky_array,Data,title,levels = np.zeros(0)):

    title_string = str(title)

    plt.figure()

    plt.subplot(111, aspect='equal')
    if len(levels):
        plt.contourf(kx_array,ky_array,Data,list(levels))
    else:
        plt.contourf(kx_array,ky_array,Data)
    plt.xlabel("kx")
    plt.ylabel("ky")
    plt.title(title_string)

    plt.colorbar()

    plt.show()
    plt.close()

def contourPlot3D_BZ(X,Y,Z,title,levels = np.zeros(0)):
    X,Y = np.meshgrid(X,Y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,cmap='viridis', edgecolor='none')
    ax.set_title('surface');

