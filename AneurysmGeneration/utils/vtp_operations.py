'''
	vtp_operations.py


'''

import vtk
import numpy
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


	NoP = model.GetNumberOfPoints()
	points = np.zeros((NoP, 3))
	for i in range(NoP):
		points[i] = model.GetPoint(i)

	if out_file_name is None:
		write_to_file(name[:-3] + '.pkl', points)
