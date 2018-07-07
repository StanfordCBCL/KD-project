'''
	vtp_operations.py


'''

import vtk
import numpy as np
from batch import write_to_file


def read_centerline_vtp(name, out_file_name=None):
	'''
		input:
			* 

		output:
			* 


	'''

	# open the polydata
	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(name)
	reader.Update()
	polydata = reader.GetOutput()


	NoP = polydata.GetNumberOfPoints()
	points = np.zeros((NoP, 3))
	for i in range(NoP):
		points[i] = polydata.GetPoint(i)

	print points.shape
	print 'bouta write these points to pkl'
	if out_file_name is None:
		write_to_file(name[:-4], points)


if __name__ == "__main__":
	print 'testing read_centerline_vtp'
	read_centerline_vtp('RCA_cl.vtp')