'''
	test_construction.py

	Generate a centerline vtp file with 100 points from (0, 0, -2) to (0, 0, 5) for testing purposes
'''

import vtk
import numpy as np

	
def gen_cyl_and_center():
	

	points_list = [[0.0, 0.0, i] for i in np.linspace(-2.0, 5.0, 1000)]


	# Create a vtkPoints object and store the points in it
	points = vtk.vtkPoints()

	for p in points_list:
		points.InsertNextPoint(p)

	centerpoly = vtk.vtkPolyData()

	centerpoly.SetPoints(points)

	writer2 = vtk.vtkXMLPolyDataWriter()
	writer2.SetFileName("test_model_centerline.vtp")
	writer2.SetInputData(centerpoly)
	writer2.Write()


if __name__ == "__main__":

	gen_cyl_and_center()