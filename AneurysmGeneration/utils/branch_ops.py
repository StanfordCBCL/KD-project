'''

	branch_ops.py


'''

import numpy as np
import vtk
from vtk.util import numpy_support as nps 


from interpolation import *
from normalization import *
from parser import *
from slice import *
from batch import *


def shift_branches(wall_model, wall_region, intersection, affected_face_displace, face_to_cap, face_to_points, point_connectivity, easing, boost=1.0):
	'''

	'''

	print 'Preparing to shift branches'
	print '------------------------------'


	# consolidate affected branch displacements
	for faceID, displace_list in affected_face_displace.iteritems():

		print '> for face id ', faceID, 'the displace_list contains ', len(displace_list), ' vecs'

		# compute the L2 norms of each displacement vector; select single displacement with the largest L2 norm
		# displace_idx = np.argmax([np.sqrt(x**2 + y**2 + z**2) for (x,y,z) in displace_list])
		#average_displace = np.max(np.array(displace_list), axis=0)

		# affected_face_displace[faceID] = displace_list[displace_idx]
		# apply the mean displacement in each dimension
		affected_face_displace[faceID] = np.mean(displace_list, axis=0)

	# apply displacement to affected branches and associated caps
	for faceID, displace in affected_face_displace.iteritems():
		cap_points = face_to_cap[faceID]
		vessel_points = face_to_points[faceID] - set(wall_region)
		branch_ids = cap_points.union(vessel_points)

		for pointID in branch_ids:
			cur_pt = wall_model.GetPoints().GetPoint(pointID)
			new_pt = [r + boost*dr for (r, dr) in zip(cur_pt, displace)]
			wall_model.GetPoints().SetPoint(pointID, new_pt)

		branch_tilt(wall_model, , , branch_ids)

		if easing:
			branch_easing(wall_model, intersection, vessel_points, point_connectivity)	



	print 'Done shifting branches'
	print '------------------------'
	return 


def branch_tilt(wall_model, inletIDs, original_positions, branch_ids):
	'''
		
		basically we're going to colelct all the displacement values
		and find some way to approximate the angle of this inlet plane 	

		input: 
			* wall_model, the underlying polydata
			* inletIDs, the set of pointIDs corresponding to the inlet of the branch that we're operating on
			* original_positions, the np array of shape (npoints, 3) that corresponds to unshifted positions
			* branch_ids, the set of pointIDs corresponding to the entire branch we're operating on 

	'''
	print 'Preparing to tilt branches'
	print '--------------------------'

	posts = extract_points(wall_model, pointIDs=inletIDs)

	# compute the normal of inlet pre-displacement
	original_positions -= np.mean(original_positions, axis=0)
	u_original,_,_ = np.linalg.svd(original_positions.T)
	normal_original = u_original[-1]

	# compute the normal of inlet post-displacement 
	posts -= np.mean(posts, axis=0)
	u_post,_,_ = np.linalg.svd(posts.T)
	normal_post = u_post[-1]

	# determine how much to tilt 
	tilt = np.diag(normal_post/normal_original)

	# tilt branch 
	branch = extract_points(wall_model, pointIDs=branch_ids)
	branch -= np.mean(posts, axis=0)
	new_branch = np.dot(branch, tilt.T)
	new_branch += np.mean(posts, axis=0)

	# save tilts
	for pointID, point in zip(branch_ids, branch):
		wall_model.GetPoints().SetPoint(pointID, point)


	print 'done tilting branches'
	print '---------------------'
	return


def branch_easing(wall_model, intersection, vessel_points, point_connectivity, num_iterations=4, aggress=.5):
	'''

		an iterative point-moving solution to awkward branch movements
	'''

	print 'Preparing to perform branch easing'
	print '----------------------------------'

	# iteratively apply a laplacian-like smoothing operation to affected points at the intersection
	for i in range(num_iterations):

		print '>>>> Easing iteration # ', i

		easing = {}

		# build a dictionary of local point connectivity for affected points
		for pointID in intersection:
			connected = point_connectivity[pointID]

			for pointID_2 in connected:
				if pointID_2 in vessel_points:

					if pointID_2 in easing.keys():
						easing[pointID_2].append(pointID)
					else: 
						easing[pointID_2] = [pointID]


		# shift points 
		for pointID, targets in easing.iteritems():
			cur_pt = wall_model.GetPoints().GetPoint(pointID)
			target_positions = [wall_model.GetPoints().GetPoint(target) for target in targets]
			pos_shifts = [[ti - ri for (ri, ti) in zip(cur_pt, target)] for target in target_positions]
			applied_shift = aggress*np.mean(np.array(pos_shifts), axis=0) 
			new_pt = [r + dr for (r, dr) in zip(cur_pt, applied_shift)]
			wall_model.GetPoints().SetPoint(pointID, new_pt)


		intersection = easing.keys()


	print 'Completed branch easing'
	print '-----------------------'



def organize_intersections(wall_region, point_to_face, cur_face):
	'''
		Adjusts data strctures to record the pointIDs at the intersection of branching vessels and initialize 
		storage for the average displacement of the affected faces. 
	'''
	print 'organizing intersections against', cur_face
	print '----------------------------------------'

	print len(wall_region)
	intersect = {}
	affected_face_set = set()

	for pt in wall_region:
		belong_faces = point_to_face[pt]
		if len(belong_faces) > 1:
			for face in belong_faces: 
				if face != cur_face:
					intersect[pt] = face
					affected_face_set.add(face)

	print 'organized intersections, the affected face set is ', affected_face_set
	print '--------------------------------'

	affected_face_displace = {faceID:[] for faceID in affected_face_set}

	return (affected_face_displace, intersect)


def centerline_shift(start_id, end_id, axial_pos, wall_model, cl_start, cl_end):
	'''
	'''
	start_pts = np.zeros((len(start_id), 3))
	end_pts = np.zeros((len(end_id), 3))

	for i, pointID in enumerate(start_id): start_pts[i] = wall_model.GetPoints().GetPoint(pointID)
	for i, pointID in enumerate(end_id): end_pts[i] = wall_model.GetPoints().GetPoint(pointID)

	centroid_start = np.mean(start_pts, axis=0)
	centroid_end = np.mean(end_pts, axis=0)

	print centroid_start
	print centroid_end

	disp_start = np.array(cl_start) - centroid_start 
	disp_end = np.array(cl_end) - centroid_end

	axial_pos = axial_pos.reshape(len(axial_pos), 1)

	adjust = (disp_end - disp_start)*axial_pos + disp_start 

	return adjust


# def smoothing(wall_model, point_connectivity, point_set, num_iterations=10, aggress=.2):
# 	'''
# 	'''

# 	for c in num_iterations:

# 		points_with_neighbors = set()

# 		pos_cache = np.zeros((len(point_set), 3))
# 		for i, pointID in enumerate(point_set): pos_cache[i] = wall_model.GetPoints().GetPoint(pointID)







