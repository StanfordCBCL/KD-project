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

	'''
	input:
		* x_interp, a vector normalized centerline coordinates that each need to be assigned an expansion coefficient 
		* centerrange, a tuple that describes the range of x coordinates that must be assigned 
		* rad_shape, an optional descriptor that allows for different shapes of aneurysms to be specified

	output: 
		* interpolated, a vector with len = len(x_interp) that contains radial expansion coefficients
		* 

	interpolated_points returns a set of radial expansion coefficients that fit a specificied shape by using scipy's 
	spline interpolation to first represent a set of points with spline, and then return the predicted fits corresponding 
	to new points. 
		i.e. we put in a shape as vector x, vector y, and then another vector my_x to get a vector my_y .

	We will be able to specify the shape of the aneurysm depending on our initial shape governed by x, y. 
	'''


	print '------------------------'
	print 'interpolating points'

	
	#initialize variable to hold the number of points provided in our predefined shape 
	piecewise_NoP = -1

	#default expansion shape; specify by specifying a set of y coords that we assume to be evenly spaced
	if rad_shape is None:
		print 'default rad_shape will be used'
		#rad_shape = [0, .1, .3, .5, .6, .6, .5, .3, .1, 0]
		rad_shape = [0, .25, .5, 1, 1, .5, .25, 0]
		piecewise_NoP = len(rad_shape)


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


