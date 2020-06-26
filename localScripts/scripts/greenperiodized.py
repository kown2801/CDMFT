import numpy as np
import json
import os
import glob
import re
clusterSize = 4
pi = np.pi

def merged(A,B,C,D):
    temp1 = np.concatenate((A,B),axis=1)
    temp2 = np.concatenate((C,D),axis=1)
    return np.concatenate((temp1,temp2),axis=0)

def split(A,size):
    temp = np.array((A[:size,:size],A[:size,size:],A[size:,:size],A[size:,size:]))
    return temp

def get_max_iter(folder):
    results_list = []
    for f in glob.glob(os.path.join(folder,"self*.dat")):
        results_list.append(int(re.search(r"self([0-9]+)\.dat",f).group(1)))
    try:
        return(max(results_list))
    except:
        print("There is no output to delete, take a look")
        exit()
        
def load_matsubara_measurements_from_simulations(dataFolder,measurement_name,converged_from,errors=False):
    max_iter = get_max_iter(dataFolder)
    i=1
    G = []
    for i in range(max_iter-converged_from+1,max_iter+1):
        G_total = np.loadtxt(os.path.join(dataFolder,measurement_name + str(i) + ".dat"))
        G.append(G_total[:,1:])
        omegas = G_total[:,0]
    #Now get the complex
    G = np.array(G,dtype=np.complex)
    G = G[:,:,::2] + 1j*G[:,:,1::2]
    return G

def get_link_matrix_from_file(folder):
    f = open(folder + "IN/LinkN.json")
    linkN_matrix = np.array(json.loads(f.read()))
    f.close()
    f = open(folder + "IN/LinkA.json")
    linkA_matrix = np.array(json.loads(f.read()))
    f.close()
    return linkA_matrix,linkN_matrix
        
    
def get_from_link(green_file_input_at_omega,linkN,linkA):
    green_up = np.zeros((clusterSize,clusterSize),dtype=np.complex)
    green_side = np.zeros((clusterSize,clusterSize),dtype=np.complex)
    dico = {"00":0,"01":1,"11":2,"pphi":3,"mphi":4,"empty":5}
    green_file_input_at_omega = np.append(green_file_input_at_omega, 0)
    for i,row in enumerate(green_up):
        for j,val in enumerate(row):
            green_up[i,j] = green_file_input_at_omega[dico[linkN[i,j]]]
    for i,row in enumerate(green_side):
        for j,val in enumerate(row):
            green_side[i,j] = green_file_input_at_omega[dico[linkA[i,j]]]
    return merged(green_up,green_side,np.conj(np.transpose(green_side)),-np.conj(green_up))

def construct_green_matrix(green_file_input_at_omega):
    linkA,linkN = get_link_matrix_from_file("Energy0/")
    green_constructed = get_from_link(green_file_input_at_omega,linkN,linkA)
    return green_constructed


def construct_full_green_matrix(green_from_file):
    G = np.zeros((green_from_file.shape[0],2*clusterSize,2*clusterSize),dtype=complex)
    for i,element in enumerate(green_from_file):
        G[i] = construct_green_matrix(element)
    return G
        


def periodize(number_of_k, green_function):
    numberOmega = green_function.shape[0]
    numberPoints = number_of_k

    positionMatrix = [[0, 0], [1, 0], [1, 1], [0, 1]]  # position of the sites 1,2,3,4 in the cluster
    sites = range(0, clusterSize)

    deltaX = np.zeros((clusterSize, clusterSize))
    deltaY = np.zeros((clusterSize, clusterSize))

    for i in sites:
        for j in sites:
            deltaX[i][j] = positionMatrix[i][0] - positionMatrix[j][0]
            deltaY[i][j] = positionMatrix[i][1] - positionMatrix[j][1]

    kx_array = np.linspace(-pi, pi, numberPoints)
    ky_array = np.linspace(-pi, pi, numberPoints)

    rangeK = range(0, numberPoints)
    K_reciprocre = [[2 * pi, 2 * pi], [2 * pi, pi], [pi, pi], [pi, 2 * pi]]

    GreenPeriodized = np.zeros((numberPoints, numberPoints, numberOmega, 2, 2), dtype="complex128")
    #print("\t\t Periodization is being made...")
    #pbar = iProgressBar(len(rangeK))

    for kx in rangeK:

        for ky in rangeK:

            for w in range(0, numberOmega):

                GreenNambu =  green_function[w]

                PartsNambu = split(GreenNambu, clusterSize)

                temp = np.exp(1j * (deltaX * kx_array[kx] + deltaY * ky_array[ky]))
                for i in range(0, len(PartsNambu)):
                    PartsNambu[i] *= temp

                GreenPeriodized[kx][ky][w][0][0] = (np.sum(PartsNambu[0])) / clusterSize
                GreenPeriodized[kx][ky][w][0][1] = (np.sum(PartsNambu[1])) / clusterSize
                GreenPeriodized[kx][ky][w][1][0] = (np.sum(PartsNambu[2])) / clusterSize
                GreenPeriodized[kx][ky][w][1][1] = (np.sum(PartsNambu[3])) / clusterSize

        #pbar.increment()

    #pbar.finish()

    return GreenPeriodized