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
from utils.parser import *
from utils.slice import *

def obtain_expansion_region(wall, centerline, included_points, start=.1, length=.1):
	'''
	input: 
		* wall vtk polydata
		* centerline as np array of xyz points
		* optional region specification

	output: 
		* list of point IDs of wall points within region
		* list of point IDs of centerline points within the region
		* dictionary of wall point ID -> closest centerline point

	Obtains the IDs of wall points and center points that are within the desired region. 
	'''

	# project wall points to closest centerline point 
	normalized_wall, normalized_center, wall_to_center = projection(wall, centerline, included_points)

	# initialize datastructures
	wall_region_id = []
	axial_pos = []
	center_region_id = []	

	# determine how many points to iterate over
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center = len(centerline)

	# find the wall points projected into the desired centerline region
	for i in included_points:
		if (normalized_wall[i] >= start) and (normalized_wall[i] <= start + length): 
			wall_region_id.append(i)
			axial_pos.append(normalized_wall[i])

	# find the center points projected within the desired region
	for i in range(NoP_center):
		if (normalized_center[i] >= start) and (normalized_center[i] <= start + length):
			center_region_id.append(i)

	return (wall_region_id, center_region_id, axial_pos, wall_to_center)




def grow_aneurysm(wall_name, centerline, included_points, start=.1, length=.1):
	'''
	input: 
		* name of wall vtp file
		* np array of centerline points 
		* optional region specification

	output:
		* writes a file to current working directory named "modified_"+wall_name

	Given an input wall and centerline, artificially grow and aneurysm at the desired region. 
	'''

	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()


	# centerreader = vtk.vtkXMLPolyDataReader()
	# centerreader.SetFileName(centerline_name)
	# centerreader.Update()
	# centerline_model = centerreader.GetOutput()

	wall_region, center_region, axial_pos, wall_to_center = obtain_expansion_region(wall_model, centerline, included_points, start, length)

	expand = interpolated_points(axial_pos, (min(axial_pos), max(axial_pos)) )

	for i, wall_id in enumerate(wall_region):

		cur_pt = wall_model.GetPoints().GetPoint(wall_id)
		normal = [r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[wall_id]) ]

		new_pt = [r + expand[i]*dn for (r,dn) in zip(cur_pt, normal)]

		wall_model.GetPoints().SetPoint(wall_id, new_pt)

	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall_model)
	new.SetFileName(wall_name[:-4] + '_modified.vtp')
	new.Write()


def main():

	#wall = "test_model.vtp"
	#centerline = "test_model_centerline.vtp"

	start = .3
	length = .3 

	# define the location of models, centerlines, metadata and specify the wall_name
	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"

	# find the centerline files within the model directory and represent them as np arrays; 
	# find the names of the centerline files (without the .pth file ending)
	# note: this matches centerline name to the np array with all the point data
	centers, names = gather_centerlines(model_dir)

	# find the face IDs assigned to the cells in the model corresponding to the centerline names in the directory
	# note: this matches centerline name against its faceID 
	corresponding_faces, face_list = parse_facenames(names, model_dir)

	# find the points corresponding to each relevant face ID
	# note: face_to_points is from faceID to list of pointID
	nonAortaPts, face_to_points = wall_isolation(face_list, model_dir=model_dir, wall_name=wall_name)

	# prepare to grow an aneurysm at a specific region along a specified centerline
	# to do this, we input the set of centerline points as an np array of [xyz] and 
	# the set of pointIDs corresponding to the right wall region
	# 
	cur_name = names[0]
	cur_face = corresponding_faces[cur_name]
	cur_center = centers[cur_name]
	cur_points = face_to_points[cur_face]
	grow_aneurysm(wall_name, cur_center, cur_points, start=start, length=length)

	# normalize coordinates for each centerline
	# norm_centers = normalized_all_centerlines(centers)

	# assign each point the normalized coordinate for its closest projected centerline point


	

if __name__ == "__main__":

	main()





