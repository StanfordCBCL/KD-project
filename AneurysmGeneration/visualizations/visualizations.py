'''
visualization.py
'''

import vtk
from vtk.util import numpy_support as nps 
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_centerlines(names, centers, resampled):
	'''
		input:
			* names, a list of the centerline names
			* centers, a dictionary of {name: list_of_points}
			* resampled, a list of resampled centerline points

		output: 

		Does a 3D scatter plot of all centerlines pre/post resampling
	'''

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for centerline in resampled:
		ax.scatter(centerline[:,0], centerline[:,1], centerline[:,2])
	for name in names:
		center = centers[name]
		ax.scatter(center[:,0], center[:,1], center[:,2])
	plt.show()


def visualize_total(name="total_cl.vtp"):

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(name)
	reader.Update()
	model = reader.GetOutput()

	NoP = model.GetNumberOfPoints()
	points = np.zeros((NoP, 3))
	for i in range(NoP):
		points[i] = model.GetPoint(i)
	print points

	ax.scatter(points[:,0], points[:, 1], points[:, 2])
	plt.show()


if __name__ == "__main__":
	VMTK_centerline()
