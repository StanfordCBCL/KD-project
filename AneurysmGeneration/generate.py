'''
	generate.py

	Given an input wall vtp and centerline vtp, distort the wall to artificially generate
	an aneurysmal region. 

'''
import numpy as np
import vtk
from vtk.util import numpy_support as nps 


from utils.interpolation import *
from utils.normalization import *


def obtain_expansion_region(wall, centerline, start=.1, length=.1):
	'''
	input: 
		* wall vtk polydata
		* centerline vtk polydata
		* optional region

	output: 
		* list of point IDs of wall points within region
		* list of point IDs of centerline points within the region
		* dictionary of wall point ID -> closest centerline point

	Obtains the IDs of wall points and center points that are within the desired region. 
	'''

	# project wall points to closest centerline point 
	normalized_wall, normalized_center, wall_to_center = projection(wall, centerline)

	# initialize datastructures
	wall_region_id = []
	axial_pos = []
	center_region_id = []	

	# determine how many points to iterate over
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center = centerline.GetNumberOfPoints()

	# find the wall points projected into the desired centerline region
	for i in range(NoP_wall):
		if (normalized_wall[i] >= start) and (normalized_wall[i] <= start + length): 
			wall_region_id.append(i)
			axial_pos.append(normalized_wall[i])

	# find the center points projected within the desired region
	for i in range(NoP_center):
		if (normalized_center[i] >= start) and (normalized_center[i] <= start + length):
			center_region_id.append(i)

	return (wall_region_id, center_region_id, axial_pos, wall_to_center)




def grow_aneurysm(wall_name, centerline_name, start=.1, length=.1):
	'''
	input: 
		* name of wall vtp file
		* name of centerline vtp file
		* optional region specification

	output:
		* writes a file to current working directory named "modified_"+wall_name

	Given an input wall and centerline, artificially grow and aneurysm at the desired region. 
	'''

	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

	centerreader = vtk.vtkXMLPolyDataReader()
	centerreader.SetFileName(centerline_name)
	centerreader.Update()
	centerline_model = centerreader.GetOutput()

	wall_region, center_region, axial_pos, wall_to_center = obtain_expansion_region(wall_model, centerline_model, start, length)

	expand = interpolated_points(axial_pos, (min(axial_pos), max(axial_pos)) )

	for i, wall_id in enumerate(wall_region):

		cur_pt = wall_model.GetPoints().GetPoint(wall_id)
		normal = [r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[wall_id]) ]

		new_pt = [r + expand[i]*dn for (r,dn) in zip(cur_pt, normal)]

		wall_model.GetPoints().SetPoint(wall_id, new_pt)

	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall_model)
	new.SetFileName("modified_"+wall_name)
	new.Write()


def main():

	#wall = "test_model.vtp"
	#centerline = "test_model_centerline.vtp"

	start = .3
	length = .3 

	# define the location of models, centerlines, metadata
	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"

	# find the centerline files within the model directory and represent them as np arrays; 
	# find the names of the centerline files (without the .pth file ending)
	centers, names = gather_centerlines(model_dir)

	# find the face IDs assigned to the cells in the model corresponding to the centerline names in the directory
	corresponding_faces = parse_facenames(names, model_dir)

	# find the points corresponding to each relevant face IDs


	# normalize coordinates for each centerline


	# assign each point the normalized coordinate for its closest projected centerline point

	
	grow_aneurysm(wall, centerline, start = start, length=length)


if __name__ == "__main__":

	main()





