'''
	interpolation.py

	Given a specified normalized centerline region bounded by 
		(a, b) s.t. a, b in (0.0, 1.0) and a < b, 
	provide an interpolated expansion coefficient for each wall point provided. 
	The shape of the expansion function can be specified. 
	Boundary conditions can also be specified. 
'''

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from pathreader import *

def interpolated_points(x_interp, centerrange, rad_shape=None, interp_type='cubic_clamped'):

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


	
	print 'interpolating points'
	print '------------------------'


	#initialize variable to hold the number of points provided in our predefined shape 
	piecewise_NoP = -1

	#initialize interpolated values
	interpolated = None

	#default expansion shape; specify by specifying a set of y coords that we assume to be evenly spaced
	if rad_shape is None:
		print 'default rad_shape will be used'
		#rad_shape = [0, .1, .3, .5, .6, .6, .5, .3, .1, 0]
		# rad_shape = [0, .75, 1.5, 2, 2, 1.5, .75, 0]
		rad_shape =[0, 1, 2, 3, 4, 4, 3, 2, 1, 0]

	else:
		print 'rad_shape has been specified:'
		start, middle, end = rad_shape
		rad_shape = [start, .5*(start + middle), middle, .5*(middle + end), end]

	# if it's not none, expect it to be a tuple of (start prescribed radius, max prescribed radius, end prescribed radius)
	# so we should maybe supplement it with some other points

	piecewise_NoP = len(rad_shape)
	print 'the rad shape to be used for aneurysm growth:'
	print rad_shape

	if interp_type is 'cubic_clamped':
		print 'cubic spline interpolation will be used with clamped bc'
		cs = interpolate.CubicSpline(np.linspace(centerrange[0], centerrange[1], piecewise_NoP), rad_shape, bc_type='clamped')
		interpolated = cs(x_interp)

	elif interp_type is 'b-spline':
		print 'b-spline interpolation will be used'
		tck = interpolate.splrep(np.linspace(centerrange[0], centerrange[1], piecewise_NoP), rad_shape)
		interpolated = interpolate.splev(x_interp, tck, der=0)

	print 'done interpolating points '
	print '--------------------------'
	return interpolated


def resample_centerline(centerline, length=None):

	print 'trying to resample centerline'
	print centerline.shape
	# spline parameters
	s=0.01 	# smoothness parameter
	k=3 	# spline order
	nest=-1 # estimate of number of knots needed (-1 = maximal)
	nSpoints_factor = 30

	centerline = np.transpose(centerline)

	# find the knot points
	tckp,u = interpolate.splprep(centerline, s=s,k=k)


	print centerline.shape

	if length is None:
		length = len(centerline[0])
	nSpoints = length*30

	print '-----   the number of spline points is ', nSpoints

	# evaluate spline, including interpolated points
	xs,ys,zs = interpolate.splev(np.linspace(0,1,nSpoints),tckp)

	print 'done resampling centerline'
	print '--------------------------'
	#print np.stack((xs, ys, zs)).shape
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
