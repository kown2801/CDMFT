#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5
from subprocess import run
import numpy as np
import argparse
import sys
import os
import glob
import re
file_dir_from_main = "DataDifferentparameters/Analysis_code/"

single_occupation_files = ["ekin.dat","pn.dat"]
single_filenames = ["Chi0.dat","Chi0Sites.dat","D.dat","DSites.dat","ekin.dat","k.dat","kSites.dat","N.dat","NSites.dat","pn.dat","sign.dat","Sz.dat","SzSites.dat"]
multiple_filenames = ["ChiFull","ChiFullSites","dgreen","green","hyb","pK","pxgreen","pxygreen","pygreen","self"]

#Beginning function to handle uncomputed occupations
def find_missing_occupation(input_dir,data_dir):
	lst = np.sort(np.loadtxt(os.path.join(data_dir,"ekin.dat"))[:,0].astype(int))
	should_be = np.sort(np.loadtxt(os.path.join(data_dir,"N.dat"))[:,0].astype(int))
	for i in should_be:
		if i not in lst:
			return i
	return None
def all_occupations_in_folder(files_dir):
	autocoherence_dir = "Autocoherence/"
	solver_dir = "ImpuritySolver/"
	input_dir = os.path.join(files_dir,"IN/")
	output_dir = os.path.join(files_dir,"OUT/")
	data_dir = os.path.join(files_dir,"DATA/")
	all_data_dir = os.path.join(files_dir,"../")
	index = find_missing_occupation(input_dir,data_dir)
	if index != None:
		run([os.path.join(files_dir,"../run_occupation.sh"), autocoherence_dir, output_dir, data_dir, "params", str(index)])
		run(["sbatch",os.path.join(file_dir_from_main,"actions.sh"),"-a","occupations","-f",os.path.basename(files_dir),"-m"])
	else:
		for file in single_occupation_files:
			data = np.loadtxt(os.path.join(data_dir,file))
			data = data[np.argsort(data[:,0])]
			np.savetxt(os.path.join(data_dir,file),data,fmt='%i %f')
		#We need to sort everything to be sure it's ok
	return 0
def order_and_remove_occupation(files_dir):
	autocoherence_dir = "Autocoherence/"
	solver_dir = "ImpuritySolver/"
	input_dir = os.path.join(files_dir,"IN/")
	output_dir = os.path.join(files_dir,"OUT/")
	data_dir = os.path.join(files_dir,"DATA/")
	all_data_dir = os.path.join(files_dir,"../")
	for file in single_filenames:
		data = np.loadtxt(os.path.join(data_dir,file))
		treated = data[np.unique(data[:,0],return_index = True)[1],:]
		try:
			np.savetxt(os.path.join(data_dir,file),treated,fmt="%i %f")
		except:
			np.savetxt(os.path.join(data_dir,file),treated,fmt="%i %f %f %f %f")
	return 0
def handle_occupation(files_dir,mode):
	if not mode:
		order_and_remove_occupation(files_dir)
	all_occupations_in_folder(files_dir)
#END

#Beginnning functions to delete the last iteration
def find_last_results(files_dir):
	results_list = []
	for f in glob.glob(os.path.join(files_dir,"OUT/params*.meas.json")):
		results_list.append(int(re.search(r"params([0-9]+)\.meas\.json",f).group(1)))
	try:
		return(max(results_list))
	except:
		print("There is no output to delete, take a look")
		exit()
def delete_safely(file):
	try:
		os.remove(file)
	except Exception as e:
		print("Erreur au fichier " + file)
def delete_last(files_dir,mode=True):
	iteration = find_last_results(files_dir)
	data_dir = os.path.join(files_dir,"DATA/")
	input_dir = os.path.join(files_dir,"IN/")
	output_dir = os.path.join(files_dir,"OUT/")
	for file in single_filenames:
		try:
			data = np.loadtxt(os.path.join(data_dir,file))
			if int(data[-1,0]) == int(iteration):
				try:
					np.savetxt(os.path.join(data_dir,file),data[:-1,:],fmt="%i %f")
				except:
						np.savetxt(os.path.join(data_dir,file),data[:-1,:],fmt="%i %f %f %f %f")
			else:
				print("Erreur au fichier " + file + " " + str(int(data[-1,0])) + " :  " + str(int(iteration)))
		except Exception as e:
			print("Erreur au fichier " + file)
	for file in multiple_filenames:
		delete_safely(os.path.join(data_dir,file + str(iteration) + ".dat"))
	delete_safely(os.path.join(input_dir,"Hyb" + str(int(iteration)+1) + ".json"))	
	delete_safely(os.path.join(input_dir,"params" + str(int(iteration)+1) + ".json"))		
	delete_safely(os.path.join(output_dir,"params" + str(iteration) + ".meas.json"))	
	return 0
#END

actions_list = ["occupations","delete_last"]
actions_functions = [handle_occupation,delete_last]

#Activate the parser to understand the input
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-v","--verbose", action="store_true",help="increase output verbosity")
parser.add_argument("-a",dest="action",help="the action of the program among : \n\t" + str(actions_list),required=True)
parser.add_argument("-f",dest="folder",help="the folder you want to act on",required=True)
parser.add_argument("-m",dest="mode",action="store_true",help="specified when computing the occupation in order to not remove doubles") 
args = parser.parse_args()
if args.action not in actions_list:
	parser.error("-a should be in " + str(actions_list))

#Get the function we want to launch
index_in_list = actions_list.index(args.action)
#Initialize the folder structure correctly
path_to_main_dir = "../../"
path_from_main_to_data = "DataDifferentparameters/"
os.chdir(os.path.join(os.path.dirname(os.path.realpath(__file__)),path_to_main_dir))
path_from_main_to_files = os.path.join(path_from_main_to_data + args.folder)
#Launch the action we want
actions_functions[index_in_list](path_from_main_to_files,args.mode)