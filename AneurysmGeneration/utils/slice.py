'''
	slice.py
'''

import numpy as np
import vtk
from parser import *

def slice_wall(wall, centerline, start=.1, length=.1):
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