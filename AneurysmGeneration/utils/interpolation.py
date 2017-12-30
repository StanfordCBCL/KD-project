import numpy as np
from scipy import interpolate


def interpolated_points(x_interp, centerrange, rad_shape=None):

	if rad_shape is None:
		rad_shape = [0, .5, 1, .5, 0]

	print np.linspace(centerrange[0], centerrange[1], 3)
	print rad_shape

	tck = interpolate.splrep(np.linspace(centerrange[0], centerrange[1], 5), rad_shape, s=0)

	interpolated = interpolate.splev(x_interp, tck, der=0)

	return interpolated
