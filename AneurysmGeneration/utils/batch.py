'''
	batch.py

	supports some operations required for batch processing of aneurysms

'''
import numpy as np
import pickle


def aggregate_options(start=.1, length=.1, rad_max=.8, easing=True, expansion_mode='absolute'):
	'''
		Combine all the options into one dict for readability and such
	'''
	return {'start':start, 'length':length, 'rad_max':rad_max, 'easing':easing, 'expansion_mode':expansion_mode}


def write_to_file(name, obj):
	'''
		Write object of specified name to a pickle file 
	'''
	path = 'pickles/' + name + '.pkl'
	file = open(path, 'wb')
	pickle.dump(obj, file)
	file.close()


def read_from_file(name):
	'''
		Return loaded object given by input name
	'''
	path = 'pickles/' + name + '.pkl'
	file = open(path, 'rb')
	new_obj = pickle.load(file)
	file.close()

	return new_obj


if __name__ == "__main__":
	print 'testing batch.py'
	print '--------------------------------------'

	my_dict = {'a':1, 'b':2, 'c':3}
	print my_dict

	write_to_file('mydict', my_dict)
	new_dict = read_from_file('mydict')
	print new_dict

	print 'done testing batch.py'
	print '--------------------------------------'