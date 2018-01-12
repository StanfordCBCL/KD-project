'''
	interpolation.py

	Given a specified normalized centerline region (a, b) s.t. a, b in (0.0, 1.0) and a < b, 
	provide an interpolated expansion coefficient for each wall point provided. 

	The shape of the expansion function can be specified. 
'''

import numpy as np
from scipy import interpolate


def interpolated_points(x_interp, centerrange, rad_shape=None):

	piecewise_NoP = -1
	if rad_shape is None:
		rad_shape = [0, .5, 1, .5, 0]
		piecewise_NoP = 5

	print np.linspace(centerrange[0], centerrange[1], 3)
	print rad_shape

	tck = interpolate.splrep(np.linspace(centerrange[0], centerrange[1], piecewise_NoP), rad_shape, s=0)

	interpolated = interpolate.splev(x_interp, tck, der=0)

	return interpolated
