'''
	interpolation.py

	Given a specified normalized centerline region (a, b) s.t. a, b in (0.0, 1.0) and a < b, 
	provide an interpolated expansion coefficient for each wall point provided. 

	The shape of the expansion function can be specified. 
'''

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def interpolated_points(x_interp, centerrange, rad_shape=None):

	piecewise_NoP = -1

	if rad_shape is None:
		#rad_shape = [0, .1, .3, .5, .6, .6, .5, .3, .1, 0]
		rad_shape = [0, .25, .5, 1, 1, .5, .25, 0]
		piecewise_NoP = len(rad_shape)

	#print np.linspace(centerrange[0], centerrange[1], 3)
	print rad_shape

	tck = interpolate.splrep(np.linspace(centerrange[0], centerrange[1], piecewise_NoP), rad_shape)

	interpolated = interpolate.splev(x_interp, tck, der=0)

	return interpolated



if __name__ == "__main__":
	print 'testing interpolation'
	print '---------------------'

	x_interp = np.linspace(.3, .6, 50)
	centerrange = (.3, .6)
	result = interpolated_points(x_interp, centerrange)
	plt.plot(x_interp, result, x_interp, [(1+a)*b for (a, b) in zip(result, np.ones(50))])
	plt.show()


