'''
	batch.py

	supports some operations required for batch processing of aneurysms

'''
import numpy as np
import pickle
import os
import vtk


def return_unstructured(path, return_reader=False): 
	"""Summary
	
	Args:
	    path (TYPE): Description
	    return_reader (bool, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(path)
	reader.Update()

	unstructured = vtk.vtkUnstructuredGrid()
	unstructured = reader.GetOutput()

	if return_reader: 
		return (unstructured, reader)

	return unstructured


def return_polydata(path, return_reader=False):
	'''
		Given a file name, open and return the data
	'''

	print 'return polydata'
	
	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(path)
	reader.Update()
	polydata = reader.GetOutput()

	if return_reader:
		return (polydata, reader)

	return polydata


def aggregate_options(start=.1, length=.1, rad_max=.8, easing=True, cur_name=None, cur_center=None, cur_face=None, expansion_mode='absolute', suffix=None):
	'''
		Combine all the options into one dict for readability and such
	'''
	return {'start':start, 'length':length, 'rad_max':rad_max, 'easing':easing, 'centerline':cur_center, 'cur_face':cur_face, 'suffix':suffix, 'expansion_mode':expansion_mode}


def write_to_file(name, obj):
	'''
		Write object of specified name to a pickle file 
	'''

	print 'writing structures to pickle'
	print '----------------------------'

	path = os.getcwd() + '/AneurysmGeneration/pickles/' + name + '.pkl'
	file = open(path, 'wb')
	pickle.dump(obj, file)
	file.close()


def read_from_file(name):
	'''
		Return loaded object given by input name
	'''
	print 'reading structures from pickle'
	print '------------------------------'

	path = os.getcwd() + '/AneurysmGeneration/pickles/' + name + '.pkl'
	file = open(path, 'rb')
	new_obj = pickle.load(file)
	file.close()

	return new_obj


def read_targets(fname, as_dict=False):
	'''
	reads a file called targets.txt that contains: 
	[vessel id] [start pos] [length] [rad_max] [suffix] 
	as a line for each artificial aneurysm that we want to generate

	[suffix] is an identifier for each aneurysm and is appended to filename after creation

	'''
	targets = []
	fname = os.getcwd() + fname

	with open(fname) as f:
		all_lines = [line.split() for line in f.readlines()]
		for line in all_lines:
			if '#' in line[0]:
				continue
			vessel = int(line[0])
			suffix = line[-1]
			start, length, rad_max = [float(val) for val in line[1:4]]
			targets.append((vessel, start, length, rad_max, suffix))

	if as_dict: 
		return {target[-1]:target[:-1] for target in targets}
	return targets


def collect_target_dictionary(side='R'): 
	"""Summary
	
	Args:
	    fname (TYPE): Description
	    side (str, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	if side == 'R': 
		return {
		'ASI2': {'p1': [2, .17, 1.3068, .3267],
				'p2': [2, .17, 1.5404, .3851], 
				'p3': [2, .17, 1.7740, .4435], 
				'p4': [2, .17, 2.0076, .5019], 
				'p5': [2, .17, 2.22412, .5603], 
				'm1': [2, .55, 1.3068 ,.3267], 
				'm2': [2, .55, 1.5404, .3851], 
				'm3': [2, .55, 1.7740, .4435], 
				'm4': [2, .55, 2.0076, .5019], 
				'm5': [2, .55, 2.22412, .5603], 
				'd1': [2, .77, 1.3068 ,.3267], 
				'd2': [2, .77, 1.5404, .3851], 
				'd3': [2, .77, 1.7740, .4435], 
				'd4': [2, .77, 2.0076, .5019], 
				'd5': [2, .77, 2.22412, .5603], 
				},
		'ASI6': {'p1': [2, .17, 3.9205, .3267], 
				'p2': [2 ,.17, 4.6213, .3851 ], 
				'p3': [2, .17, 5.3221, .4435], 
				'p4': [2, .17, 6.0229, .5019], 
				'p5': [2, .17, 6.7236, .5603], 
				},
		'ASI4': { 'p1': [2, .17, 2.6137, .3267], 
				'p2': [2, .17, 3.0809, .3851], 
				'p3': [2, .17, 3.5481, .4435], 
				'p4': [2, .17, 4.0152, .5019], 
				'p5': [2, .17, 4.4824, .5603],
				}
		}
	elif side == 'L': 
		return {
		'ASI2': {'lad1': [0, .23, 1.1846, .29615], 
				'lad2': [0, .23, 1.3801, .34502], 
				'lad3': [0, .23, 1.5755, .39388],
				'lad4': [0, .23, 1.7710, .44274], 
				'lad5': [0, .23, 1.9664, .49161],
				},
		'ASI6': {'lad1': [0, .235, 3.5539, .29615], 
				'lad2': [0, .235, 4.1402, .34502], 
				'lad3': [0, .235, 4.7266, .39388], 
				'lad4': [0, .235, 5.3129, .44274], 
				'lad5': [0, .235, 5.8993, .49161],
				},
		'ASI4': {'lad1': [0, .235, 2.3692, .29615], 
				'lad2': [0, .235, 2.7601, .34502],
				'lad3': [0, .235, 3.1510, .39388 ], 
				'lad4': [0, .235, 3.5419, .44274], 
				'lad5': [0, .235, 3.9328, .49161]
				},
		}
	


def batch_targets(names, corresponding_faces, resampled, targets_name, batch_status, easing):
	'''

	'''

	targets = read_targets(targets_name)
	
	if not batch_status:
		targets = [targets[0]]
	
	options = []
	for target in targets:
		(vessel, start, length, rad_max, suffix) = target
		options.append(aggregate_options(
			start=start,
			length=length,
			rad_max = rad_max,
			easing = easing,
			cur_name = names[vessel],
			cur_face = corresponding_faces[names[vessel]],
			cur_center = resampled[vessel], 
			suffix = suffix))

	return options

	

if __name__ == "__main__":
	print 'testing batch.py'
	print '--------------------------------------'

	# my_dict = {'a':1, 'b':2, 'c':3}
	# print my_dict

	# write_to_file('mydict', my_dict)
	# new_dict = read_from_file('mydict')
	# print new_dict

	print read_targets()

	print 'done testing batch.py'
	print '--------------------------------------'