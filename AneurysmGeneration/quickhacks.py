'''
	quickhacks.py

'''
import vtk
from vtk.util import numpy_support as nps 
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

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


	# dispaly a little bit first 
	
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')
	# ax.scatter(centerline[:,0], centerline[:,1], centerline[:,2])
	# ax.scatter(centerline[600:,0], centerline[600:,1], centerline[600:,2])
	# plt.show()

	centerline = centerline[600:]


	#diff = np.diff(centerline, axis=0)
	#print diff
	#okay = np.where(np.abs(diff[:,0]) + np.abs(diff[:,1]) + np.abs(diff[:,2]) > 0)
	#print okay
	#centerline = np.r_[np.reshape(centerline[0], (1, 3)), centerline[okay], np.reshape(centerline[-1,:], (1, 3))]

	centerline = np.unique(centerline, axis=0)
	# interpolate so that we get a smoother point matching 
	centerline = resample_centerline(centerline)

	write_to_file(name, centerline)

if __name__ == "__main__":
	VMTK_centerline()
