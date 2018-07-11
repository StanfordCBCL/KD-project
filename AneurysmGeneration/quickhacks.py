'''
	quickhacks.py

'''
import vtk
from vtk.util import numpy_support as nps 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from utils.batch import read_from_file, write_to_file
from utils.interpolation import resample_centerline

def VMTK_centerline(name='RCA_cl'):
	'''
		input: 
			* name of centerline to visualize and potentially clip

		output:
			* 
	'''

	# read in the points extracted from the vtk output by vmtk script
	centerline = read_from_file(name)


	# display a little bit first 
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot3D(centerline[:,0], centerline[:,1], centerline[:,2])
	ax.scatter(centerline[600:,0], centerline[600:,1], centerline[600:,2])
	plt.show()

	centerline = centerline[600:]

	# interpolate so that we get a smoother point matching 
	_, CL_indices = np.unique(centerline, return_index=True, axis=0)
	centerline = centerline[sorted(CL_indices)]
	centerline = resample_centerline(centerline)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot3D(centerline[:,0], centerline[:,1], centerline[:,2])
	plt.show()

	write_to_file(name, centerline)

if __name__ == "__main__":
	VMTK_centerline()
