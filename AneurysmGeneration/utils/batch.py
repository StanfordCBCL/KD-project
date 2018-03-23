'''
	batch.py



'''
import numpy as np


def aggregate_options(start=.1, length=.1, rad_max=.8, easing=True, expansion_mode='absolute'):
	'''
		Combine all the options into one dict for readability and such
	'''
	return {'start':start, 'length':length, 'rad_max':rad_max, 'easing':easing, 'expansion_mode':expansion_mode}

