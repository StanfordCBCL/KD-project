'''
	normalization.py

	Provide utility functinos to project points from the wall polydata to the nearest point in centerline polydata. 

	Normalization and projections adapted from Justin's scripts. 
'''

import numpy as np
import vtk 


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

def acquire_radii(pointIDs, wall_model, wall_to_center):
	'''
		Collect the radial component distances from wall points to corresponding centerline points

		input: 
			* a list, of pointIDs used to identify points on the wlal
			* a polydata 
			* a dictionary that returns the closest centerline point for an input pointID

		output: 
			* a list of L2 norms

	'''
	radii = []
	for pointID in pointIDs:
		cur_pt = wall_model.GetPoints().GetPoint(pointID)
		radial_component = np.linalg.norm([r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[pointID]) ])
		radii.append(radial_component)

	return radii

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

	print 'normalizing the centerline'
	print '--------------------------'
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

def compute_reference_norm(centerline):

	NoP = centerline.shape[0]
	p0 = np.roll(centerline, shift = -1, axis= 1)
	p1 = np.roll(centerline, shift = 0, axis = 1)
	p2 = np.roll(centerline, shift = 1, axis = 1)

	t21 = p2 - p1
	t21_normed = np.divide(t21, np.linalg.norm(t21, axis=1).reshape(NoP, 1))

	t10 = p1 - p0
	t10_normed = np.divide(t10, np.linalg.norm(t10, axis=1).reshape(NoP, 1))

	n1 = t21_normed - t10_normed
	return np.divide(n1, np.linalg.norm(n1, axis=1).reshape(NoP, 1))

def projection(wall, centerline, included_points):
	'''
		input: 
			* wall polydata 
			* centerline points as np array of shape (NoP, 3) 
			* included_points, the list of point IDs for the vessel wall we are considering

		output:
			* list of normalized centerpoint position for each wall point 
			* list of normalized centerpoint positions 
			* dictionary of correspondences between wall point index -> closest centerpoint

		For each wall point, go through all the centerline points and find the closest one. 

		Record the closest normalized centerline distance and centerline point's coordinates. 

	'''

	print 'projecting wall points onto the centerline'
	print '------------------------------------------'
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center, normalized_center, centerline_length = normalized_centerline_pth(centerline)
	reference_norms = compute_reference_norm(normalized_center)

	print '----     centerline length:   -------'
	print '----     ', centerline_length, '     -----'
	print '-------------------------------------'
	normalized_wall = np.zeros((NoP_wall))

	wall_to_center = {}
	wall_to_norm = {}

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

	return (normalized_wall, normalized_center, wall_to_center,  centerline_length)


