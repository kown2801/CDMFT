#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/python/3.5.4/bin/python3.5
import glob, os
import sys

def delete_safely(file):
    try:
        os.remove(file)
    except:
        pass

def delete_useless(folder_name):
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    os.chdir(os.path.join("../",folder_name))
    for f in glob.glob("slurm*"):
        os.remove(f)
    delete_safely("occupation-log.out")
    delete_safely("logfile")
    delete_safely("run.sh")
    os.chdir("OUT/")
    for f in glob.glob("config_*.json"):
        os.remove(f)


if __name__ == "__main__":
    if(len(sys.argv) >= 2):
        folder = sys.argv[1]
    delete_useless(folder)