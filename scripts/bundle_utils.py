from pathlib import Path
import json
#Two classes and a two functions for json pretty printing
#They work by replacing all the lists in the dictionary by MarkedLists and specifying a way of serializing the MarkedList object

def delete_safely(file):
	try:
		file.unlink()
	except Exception as e:
		print("Erreur au fichier " + file)

class MarkedList:
	_list = None
	def __init__(self, l):
		self._list = l
class CustomJSONEncoder(json.JSONEncoder):
	def default(self, o):
		if isinstance(o, MarkedList):
			return "##<{}>##".format(o._list)
def lists_to_MarkedLists(d):
	for k, v in d.items():
		if isinstance(v, dict):
			lists_to_MarkedLists(v)
		elif isinstance(v, list):
			d[k] = MarkedList(v)
	return d
            
def dump_to_pretty_json(dico):
	dico = lists_to_MarkedLists(dico)
	temp_string = json.dumps(dico,sort_keys=True,indent='\t',cls=CustomJSONEncoder)
	return temp_string.replace('"##<', "").replace('>##"', "")

def write_json_to_file(data,filename):
	to_write = dump_to_pretty_json(data)
	filename = Path(filename)
	with filename.open('w+') as file:
		file.write(to_write)
		
#This function bundles all files of the type `file + "*.json"` in one big file, for * an integer less than max_iteration
def bundle_up_to(filename_start,filename_end,max_iteration):
	#First we need to load or create the bundled file
	bundled_filename = Path(str(filename_start) + filename_end)
	bundled_data = {}
	if bundled_filename.is_file():
		with bundled_filename.open() as file:
			bundled_data = json.load(file)
	#Then we check every file that begins with `filename_start`		
	for i in range(max_iteration):
		filename = Path(str(filename_start) + str(i) + filename_end)
		if filename.is_file(): #If the file fits the criteria, we bundle it and remove it 
			with filename.open() as file:
				this_data = json.load(file)
				bundled_data[str(i)] = this_data
			filename.unlink()
	write_json_to_file(bundled_data,bundled_filename)
		
#This function debundles a specific iteration from a specific file
def debundle(filename_start,filename_end,iteration):
	bundled_filename = Path(str(filename_start) + filename_end)
	if bundled_filename.is_file():
		bundled_data = {}
		with bundled_filename.open() as file:
			bundled_data = json.load(file)
		if str(iteration) in bundled_data:
			filename = Path(str(filename_start) + str(iteration) + filename_end)
			write_json_to_file(bundled_data[str(iteration)],filename)
			return True
	return False
#This function accesses the content in a bundle at an iteration
def get_in_bundle(filename_start,filename_end,iteration):
	bundled_filename = Path(str(filename_start) + filename_end)
	if bundled_filename.is_file():
		with bundled_filename.open() as file:
			bundled_data = json.load(file)
			if str(iteration) in bundled_data:
				return bundled_data[str(iteration)]
	return None
	
#This function allows to retrieve a file or its bundled version
def get_file_content(filename_start,filename_end,iteration):
	filename = Path(str(filename_start) + str(iteration) + filename_end)
	if filename.is_file():
		with filename.open() as file:
			return json.load(file)
	else:
		return get_in_bundle(filename_start,filename_end,iteration)

#This function allows to delete a file in the bundle
def delete_in_bundle(filename_start,filename_end,iteration):
	bundled_filename = Path(str(filename_start) + filename_end)
	if bundled_filename.is_file():
		bundled_data = {}
		with bundled_filename.open() as file:
			bundled_data = json.load(file)
		bundled_data.pop(str(iteration),None)
		write_json_to_file(bundled_data,bundled_filename)
		
#This function allows to delete a file and make the last bundled file available
def delete_and_debundle(filename_start,filename_end,iteration):
	#First we delete the file
	filename = Path(str(filename_start) + str(iteration) + filename_end)
	delete_safely(filename)	
	#Then we debundle the last one
	if debundle(filename_start,filename_end,iteration-1):
		#And delete in the bundle
		delete_in_bundle(filename_start,filename_end,iteration-1)
		
def get_iterations_in_bundle(filename_start,filename_end):
	bundled_filename = Path(str(filename_start) + filename_end)
	if bundled_filename.is_file():
		with bundled_filename.open() as file:
			bundled_data = json.load(file)
			return bundled_data.keys()
	return None
	
def get_bundle(filename_start,filename_end):
	bundled_filename = Path(str(filename_start) + filename_end)
	if bundled_filename.is_file():
		with bundled_filename.open() as file:
			return json.load(file)
	return None
	
def get_last_iterations_from_bundle(filename_start,filename_end,converged_from=None,max_iter=None):
	bundle_data = get_bundle(filename_start,filename_end)
	if converged_from == None:
		return bundle_data
	else:
		if max_iter == None:
			max_iter = max(map(int,bundle_data.keys()))
		all_data = [bundle_data[str(i)] for i in range(max_iter - converged_from+1,max_iter+1)]
		return all_data