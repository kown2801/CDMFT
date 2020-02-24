#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.8.0/bin/python3.8
import os
from subprocess import call
import json
import shutil
import sys
import subprocess

def compute_all_missing_occupations():
	path_to_main_dir = "../.."
	backup_path = "./BACKUP_START"
	supra_path = os.path.join(os.getcwd(), path_to_main_dir)
	files_path = supra_path + "/DataDifferentparameters/"
	for root, dirs, files in os.walk(files_path):
		for dir_name in dirs:
			call(["sbatch",os.path.join(files_path,"run_all_occupations.sh"),files_path,supra_path,os.path.join(files_path,dir_name)])
		break
if __name__ == "__main__":
	compute_all_missing_occupations()