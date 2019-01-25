'''
	normalization.py

	Provide utility functions to project points from the wall polydata to the nearest point in centerline polydata. 

	Normalization and projections adapted from Justin's scripts. 
'''
import sys
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


def minimize_distances(points_results, points_source):

	distances = np.sqrt(((points_results - points_source[:, np.newaxis])**2).sum(axis=2))
	return (np.argmin(distances, axis=0), np.amin(distances, axis=0))


def compute_reference_norm(centerline):
	"""Summary
	
	Args:
	    centerline (ndarray): points in the centerline
	
	Returns:
	    TYPE: Description
	"""
	print 'computing reference norms'
	print '-------------------------'

	NoP = centerline.shape[0]
	p0 = np.roll(centerline, shift = -1, axis= 1)
	p1 = np.roll(centerline, shift = 0, axis = 1)
	p2 = np.roll(centerline, shift = 1, axis = 1)

	t21 = p2 - p1
	t21_normed = np.divide(t21, np.linalg.norm(t21, axis=1).reshape(NoP, 1))

	t10 = p1 - p0
	t10_normed = np.divide(t10, np.linalg.norm(t10, axis=1).reshape(NoP, 1))

	n1 = t21_normed - t10_normed
	return (np.divide(n1, np.linalg.norm(n1, axis=1).reshape(NoP, 1)), t21)


def compute_theta(r, n, t):
	'''	

	Args:
	    r (ndarray): a vector from centerline point to wall point
	    n (ndarray): the centerline's normal vector at the centerline point corresponding to the wall point
	    t (ndarray): the tangent vector corresponding to that centerline point
	
	Returns:
	    theta (ndarray): the angle between the two vectors in [-pi, pi]
	'''
	

	# print 'computing theta'
	# print '----------------'

	# precompute some stuff
	mag_r = np.linalg.norm(r, axis=1)
	mag_n = np.linalg.norm(n, axis=1)

	# compute cos
	cos_rn = (r*n).sum(axis=1)
	cos_rn = np.divide(cos_rn, mag_r)
	cos_rn = np.divide(cos_rn, mag_n)

	# compute sin 
	cx = np.cross(r, n)
	sin_rn = np.sign((cx*t).sum(axis=1))*np.linalg.norm(cx, axis=1)/mag_r/mag_n 

	# returning arctan 
	arctan_rn = np.arctan2(sin_rn, cos_rn)

	return arctan_rn


def projection(NoP, centerline, wall_points, included_ids):
	'''
		input: 
			* wall polydata 
			* centerline, the centerline points as np array of shape (NoP_center, 3) 
			* wall_points, the points of the entire model as shape (NoP_wall, 3)
			* included_ids, the set of point IDs for the vessel wall we are considering

		output:
			* transformed_wall_ref, a np array of shape (NoP, 2) representing (axial pos, theta) for each point on the wall
			* normalized_center, a list of normalized centerpoint positions 
			* wall_to_center, dictionary of correspondences between wall point index -> closest centerpoint
			* centerline_length, the length of the centerline

		For each wall point, go through all the centerline points and find the closest one. 

		Record the closest normalized centerline distance and centerline point's coordinates. 

	'''

	print 'projecting wall points onto the centerline'
	print '------------------------------------------'

	NoP_center, normalized_center, centerline_length = normalized_centerline_pth(centerline)
	

	print '----     centerline length:   -------'
	print '----     ', centerline_length, '     -----'
	print '-------------------------------------'

	# compute the angular reference positions
	reference_norms, reference_tangents = compute_reference_norm(centerline)

	# initialize wall axial pos
	transformed_wall_ref = np.zeros((NoP, 2))

	wall_to_center = {}
	min_dists = np.zeros((NoP, 1))

	chunk_sz = 500
	n_chunks = len(included_ids)//chunk_sz

	for c in range(n_chunks):
		print '> chunk ', c, '/', n_chunks -1, '		\r',
		sys.stdout.flush()
		id_chunk = included_ids[c*chunk_sz:(c+1)*chunk_sz]
		if c == n_chunks - 1: id_chunk = included_ids[c*chunk_sz:]

		center_indices, dists = minimize_distances(wall_points[id_chunk], centerline)
		min_dists[id_chunk] = dists.reshape(len(id_chunk), 1)

		r = wall_points[id_chunk] - centerline[center_indices]
		n = reference_norms[center_indices]
		t = reference_tangents[center_indices]
		transformed_wall_ref[id_chunk,0] = normalized_center[center_indices]
		transformed_wall_ref[id_chunk,1] = compute_theta(r, n, t)

		for pointID, center_idx in zip(id_chunk, center_indices): 
			wall_to_center[pointID] = centerline[center_idx]

	print ""

	return (transformed_wall_ref, normalized_center, wall_to_center, min_dists, centerline_length)


def compute_curvature(r, rp, rpp): 
    """Compute the local curvature as given in https://en.wikipedia.org/wiki/Curvature#Local_expressions 
    
    Args:
        r (TYPE): Description
        rp (TYPE): Description
        rpp (TYPE): Description
    """

    curv = 'lol'


    return curv

if __name__ == "__main__": 

	print 'testing normalization with curv computation'

	c_path = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/wall_lca1.pth"

	centerline = np.transpose(read_centerline(c_path))

	print centerline
	lol = resample_centerline(centerline)



