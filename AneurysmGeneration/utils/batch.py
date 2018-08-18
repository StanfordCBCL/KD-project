'''
	batch.py

	supports some operations required for batch processing of aneurysms

'''
import numpy as np
import pickle
import os


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


def read_targets(fname='/AneurysmGeneration/targets.txt', as_dict=False):
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


def batch_targets(names, corresponding_faces, resampled, batch_status, easing):
	'''

	'''

	targets = read_targets()
	
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