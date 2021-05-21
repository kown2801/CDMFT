import os
import json
#Two classes and a two functions for json pretty printing
#They work by replacing all the lists in the dictionary by MarkedLists and specifying a way of serializing the MarkedList object

def delete_safely(file):
	try:
		os.remove(file)
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
	with open(filename,'w+') as file:
		file.write(to_write)
		
#This function bundles all files of the type `file + "*.json"` in one big file, for * an integer less than max_iteration
def bundle_up_to(filename_start,filename_end,max_iteration):
	#First we need to load or create the bundled file
	bundled_filename = filename_start + filename_end
	bundled_data = {}
	if os.path.isfile(bundled_filename):
		with open(bundled_filename) as file:
			bundled_data = json.load(file)
	#Then we check every file that begins with `filename_start`		
	for i in range(max_iteration):
		filename = filename_start + str(i) + filename_end
		if os.path.isfile(filename): #If the file fits the criteria, we bundle it and remove it 
			with open(filename) as file:
				this_data = json.load(file)
				bundled_data[str(i)] = this_data
			os.remove(filename)
	write_json_to_file(bundled_data,bundled_filename)
		
#This function debundles a specific iteration from a specific file
def debundle(filename_start,filename_end,iteration):
	bundled_filename = filename_start + filename_end
	if os.path.isfile(bundled_filename):
		bundled_data = {}
		with open(bundled_filename) as file:
			bundled_data = json.load(file)
		if str(iteration) in bundled_data:
			filename = filename_start + str(iteration) + filename_end
			write_json_to_file(bundled_data[str(iteration)],filename)
			return True
	return False
#This function accesses the content in a bundle at an iteration
def get_in_bundle(filename_start,filename_end,iteration):
	bundled_filename = filename_start + filename_end
	if os.path.isfile(bundled_filename):
		with open(bundled_filename) as file:
			bundled_data = json.load(file)
			if str(iteration) in bundled_data:
				return bundled_data[str(iteration)]
	return None
	
#This function allows to retrieve a file or its bundled version
def get_file_content(filename_start,filename_end,iteration):
	filename = filename_start + str(iteration) + filename_end
	if os.path.isfile(filename):
		with open(filename) as file:
			return json.load(file)
	else:
		return Bu.get_in_bundle(filename_start,filename_end,iteration)

#This function allows to delete a file in the bundle
def delete_in_bundle(filename_start,filename_end,iteration):
	bundled_filename = filename_start + filename_end
	if os.path.isfile(bundled_filename):
		bundled_data = {}
		with open(bundled_filename) as file:
			bundled_data = json.load(file)
		bundled_data.pop(str(iteration),None)
		write_json_to_file(bundled_data,bundled_filename)
		
#This function allows to delete a file and make the last bundled file available
def delete_and_debundle(filename_start,filename_end,iteration):
	#First we delete the file
	filename = filename_start + str(iteration) + filename_end
	delete_safely(filename)	
	#Then we debundle the last one
	if debundle(filename_start,filename_end,iteration-1):
		#And delete in the bundle
		delete_in_bundle(filename_start,filename_end,iteration-1)