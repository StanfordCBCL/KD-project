'''
	normalization.py

	Provide utility functinos to project points from the wall polydata to the nearest point in centerline polydata. 

	Normalization and projections adapted from Justin's scripts. 
'''

import numpy as np
import vtk
#from vtk.util import numpy_support as nps 


def normalized_centerline(centerline_model):
	'''
	input: 
		* centerline polydata 

	output:
		* number of points in the centerline, 
		* list of normalized position (0.0 to 1.0) for each point
		* total centerline length

	Iterates through the centerline points twice; assumes sum of linear distance between points is reflective
	of total length of the centerline. 
	'''

	centerline_length = 0.0
	NoP = centerline_model.GetNumberOfPoints()
	normalized = {0.0: 0.0}

 	for i in range(1, NoP):
		pt = centerline_model.GetPoints().GetPoint(i)
		pt_prev = centerline_model.GetPoints().GetPoint(i-1)
		d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, pt_prev)
		d_temp = np.sqrt(d_temp)
		centerline_length = centerline_length + d_temp

	cur_length = 0.0
	for i in range(1, NoP):
		pt = centerline_model.GetPoints().GetPoint(i)
		pt_prev = centerline_model.GetPoints().GetPoint(i-1)
		d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, pt_prev)
		d_temp = np.sqrt(d_temp)
		cur_length = cur_length + d_temp
		normalized[i] = cur_length/centerline_length


	return (NoP, normalized, centerline_length)



def projection(wall, centerline):
	'''
		input: 
			* wall polydata 
			* centerline polydata

		output:
			* list of normalized centerpoint position for each wall point 
			* list of normalized centerpoint positions 
			* dictionary of correspondences between wall point index -> closest centerpoint

		For each wall point, go through all the centerline points and find the closest one. 

		Record the closest normalized centerline distance and centerline point's coordinates. 

	'''
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center, normalized_center, centerline_length = normalized_centerline(centerline)

	normalized_wall = np.zeros((NoP_wall))

	wall_to_center = {}

	for i in range(NoP_wall):
		wall_pt = wall.GetPoints().GetPoint(i)

		min_dist = float('inf')
		min_idx = -1
		for k in range(NoP_center):
			center_pt = centerline.GetPoints().GetPoint(k)
			cur_dist = vtk.vtkMath.Distance2BetweenPoints(wall_pt, center_pt)
			if cur_dist < min_dist:
				min_dist = cur_dist
				min_idx = k

		normalized_wall[i] = normalized_center[min_idx]

		wall_to_center[i] = centerline.GetPoints().GetPoint(min_idx)

	return (normalized_wall, normalized_center, wall_to_center)


