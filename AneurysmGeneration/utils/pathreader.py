'''
	pathreader.py

	Provides functionality for parsing point data out of an XML .pth file into an ndarray of shape (n,3). 
'''

import xml.etree.ElementTree as ET 
import re 
import numpy as np
import vtk
#from parser import *
from interpolation import *
from normalization import *
from vtk.util import numpy_support as nps 


def read_centerline(path_name):
	'''
	input:
		* specify the path leading to the .pth 

	output: 
		* np array of points in shape (NoP, 3)
	'''

	# read in the .pth into a string buffer
	with open(path_name) as f:
		xml = f.read()

	# adjust the structure of the XML data to have a single root node (single top layer tag), as required for ElementTree
	root = ET.fromstring(re.sub(r"(<\?xml[^>]+\?>)", r"\1<root>", xml) + "</root>")

	# access the structure containing all path points
	path_points = root[1][0][0][1]

	# iterate through path point dictionaries and place coordinates into list
	point_list = []
	for point in path_points:
		point_coords = point[0].attrib
		xyz = [float(pos) for pos in [point_coords['x'], point_coords['y'], point_coords['z'] ]]
		point_list.append(xyz)


	# return list as np array
	return np.array(point_list)



if __name__ == "__main__":
	print "testing pathreader.py"
	print "---------------------"

	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	centers, names = gather_centerlines(model_dir)
	resampled = [resample_centerline(centers[name]) for name in names]

	print compute_reference_norm(resampled[0])

	'''
	for c, name in enumerate(names):
		cl = resampled[c]
		print cl.shape
		points = vtk.vtkPoints()

		for idx in range(cl.shape[0]):
			points.InsertNextPoint(cl[idx, :])

		#points.SetScalars(np.arange(cl.shape[0]))
		polyData = vtk.vtkPolyData()
		polyData.SetPoints(points)
		polyData.GetPointData().SetScalars(nps.numpy_to_vtk(np.arange(cl.shape[0])))

		new = vtk.vtkXMLPolyDataWriter()
		new.SetInputData(polyData)
		new.SetFileName('exist_' + name + '_cl' + '.vtp')
		new.Write()

	'''
	#read_centerline(file_loc + file_name)


