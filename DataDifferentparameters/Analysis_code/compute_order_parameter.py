#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5
import os
import math
import json
import sys
import numpy as np
import re

def compute(folder_name):
    os.chdir("../")
    filename = "green"
    save_filename = folder_name + "/order_graph"
    Composante = 3
    graph = np.zeros(0)
    G = []
    offset = 1
    try:
        graph = np.load(save_filename + ".npy")
        offset+=len(graph)
    except:
        pass
    i=offset
    try:
        while(True):
            G.append([])
            f=open(folder_name + "/DATA/" + filename + str(i) + ".dat", "r")
            fl = f.readlines()
            for num, line in enumerate(fl):
                valueList = re.split(" |\n",line)[1:-1]
                G[i-offset].append(valueList)
            i+=1
            f.close()
    except FileNotFoundError as e:
        print("Fin de la lecture des fichiers, le dernier était le numéro: " + str(i-1))
        if i-1 == 0:
            print(e)
    except Exception as e:
        print(str(e))
    if len(G[0]) != 0:
        G = np.array(G[:-1]).astype(float)
        print(G.shape)
        #Ici, G[it][iwn][site]
        G = np.transpose(G,(1,2,0)) #Ici, G[iwn][site][it], mieux pour voir quelquechose
        with open(folder_name + "/IN/params1.json") as f:
            beta = json.load(f)["beta"]
        G = np.sum(G,axis=0)*2/beta
        graph = np.concatenate((graph,(G[2*(Composante+1)] - G[2*Composante])/2))
        np.save(save_filename,graph)


if __name__ == "__main__":
    if(len(sys.argv) >= 2):
        folder = sys.argv[1]
    compute(folder)