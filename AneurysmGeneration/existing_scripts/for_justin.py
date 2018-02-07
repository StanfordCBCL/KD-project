'''
	pathreader.py

	Provides functionality for parsing point data out of an XML .pth file into an ndarray of shape (n,3). 

	I also stuck in the functions for normalizing the centerline and mapping each vessel wall point to nearest centerline point 
	which were adapted from Justin's scripts. 
'''

import xml.etree.ElementTree as ET 
import re 
import numpy as np
import vtk

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


def normalized_centerline_pth(center):
	'''
	input:
		* np array of shape (NoP, 3)

	output: 
		* NoP
		* np array of length NoP, containing normalized coordinate for each
		* total centerline length

	Assigns each centerline point a total length-normalized position, holding assigned coordinate
	in form of np array with shape (NoP,). 
	
	'''

	centerline_length = 0.0
	NoP = len(center)
	normalized = np.zeros(NoP)

	for i in range(1, NoP):
		pt = center[i]
		prev_pt = center[i-1]
		d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, prev_pt)
		d_temp = np.sqrt(d_temp)
		centerline_length += d_temp
		normalized[i] = centerline_length

	normalized /= centerline_length

	return (NoP, normalized, centerline_length)


def projection(wall, centerline, included_points):
	'''
		input: 
			* wall polydata 
			* centerline points as np array of shape (NoP, 3) 

		output:
			* list of normalized centerpoint position for each wall point 
			* list of normalized centerpoint positions 
			* dictionary of correspondences between wall point index -> closest centerpoint

		For each wall point, go through all the centerline points and find the closest one. 

		Record the closest normalized centerline distance and centerline point's coordinates. 

	'''
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center, normalized_center, centerline_length = normalized_centerline_pth(centerline)

	normalized_wall = np.zeros((NoP_wall))

	wall_to_center = {}

	for i in included_points:
		wall_pt = wall.GetPoints().GetPoint(i)

		min_dist = float('inf')
		min_idx = -1
		for k in range(NoP_center):
			center_pt = centerline[k]
			cur_dist = vtk.vtkMath.Distance2BetweenPoints(wall_pt, center_pt)
			if cur_dist < min_dist:
				min_dist = cur_dist
				min_idx = k

		normalized_wall[i] = normalized_center[min_idx]

		wall_to_center[i] = centerline[min_idx]

	return (normalized_wall, normalized_center, wall_to_center)


if __name__ == "__main__":
	print "testing pathreader.py"
	print "---------------------"

	file_loc = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	file_name = "lad.pth"

	read_centerline(file_loc + file_name)


