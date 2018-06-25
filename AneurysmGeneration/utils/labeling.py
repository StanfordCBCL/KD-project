'''

	labeling.py

	Provides functions to support visualization of point data and writing this to the output vtk
'''

import numpy as np
from vtk.util import numpy_support as nps
import vtk


'''
	interpolation.py

	Given a specified normalized centerline region (a, b) s.t. a, b in (0.0, 1.0) and a < b, 
	provide an interpolated expansion coefficient for each wall point provided. 

	The shape of the expansion function can be specified. 
	Boundary conditions can also be specified. 
'''

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from pathreader import *

def interpolated_points(x_interp, centerrange, rad_shape='expansion_coeff', rad_bounds = None, interp_type='cubic_clamped'):

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


	print centerrange
	#initialize variable to hold the number of points provided in our predefined shape 
	piecewise_NoP = -1

	# initialize the interpolated
	interpolated = None

	#default expansion shape; specify by specifying a set of y coords that we assume to be evenly spaced
	if rad_shape == 'expansion coefficient':
		print 'default rad_shape will be used'
		#rad_shape = [0, .1, .3, .5, .6, .6, .5, .3, .1, 0]
		# rad_shape = [0, .75, 1.5, 2, 2, 1.5, .75, 0]
		rad_shape =[0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
		piecewise_NoP = len(rad_shape)

	# elif rad_shape == 'specify max':
	# 	print 'specifying maximum absolute radius'
	# 	if rad_bounds is None:
	# 		print 'Please specify the initial and final radii'

	print rad_shape

	if interp_type == 'cubic_clamped':
		print 'cubic spline interpolation will be used with clamped bc'
		cs = interpolate.CubicSpline(np.linspace(centerrange[0], centerrange[1], piecewise_NoP), rad_shape, bc_type='clamped')
		interpolated = cs(x_interp)

	elif interp_type == 'b-spline':
	
		tck = interpolate.splrep(np.linspace(centerrange[0], centerrange[1], piecewise_NoP), rad_shape)
		interpolated = interpolate.splev(x_interp, tck, der=0)

	return interpolated
	# print 'oh no something went wrong'
	# return None


def resample_centerline(centerline, length=None):

	print 'trying to resample centerline'
	print centerline.shape
	# spline parameters
	s=0.01 	# smoothness parameter
	k=3 	# spline order
	nest=-1 # estimate of number of knots needed (-1 = maximal)

	centerline = np.transpose(centerline)

	# find the knot points
	tckp,u = interpolate.splprep(centerline, s=s,k=k)


	print centerline.shape

	if length is None:
		length = len(centerline[0])
	nSpoints = length*20

	# evaluate spline, including interpolated points
	xs,ys,zs = interpolate.splev(np.linspace(0,1,nSpoints),tckp)

	print np.stack((xs, ys, zs)).shape
	return np.transpose(np.stack((xs, ys, zs)))


if __name__ == "__main__":

	# print 'testing interpolation'
	# print '---------------------'

	# x_interp = np.linspace(.3, .6, 50)
	# centerrange = (.3, .6)
	# result = interpolated_points(x_interp, centerrange, interp_type='clamped')
	# plt.plot(x_interp, result, x_interp, [(1+a)*b for (a, b) in zip(result, np.ones(50))])
	# plt.show()



	print 'testing resample_centerline'
	print '---------------------'

	c_path = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/wall_lca1.pth"

	centerline = np.transpose(read_centerline(c_path))

	print centerline
	lol = resample_centerline(centerline)
	
	print lol
	





