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


def shift_branches(wall_model, wall_region, intersection, originals, affected_face_displace, face_to_cap, face_to_points, point_connectivity, easing, boost=1.0):
	'''
		inputs: 
			* wall_model, 				polydata of whole model
			* wall_region, 				list of pointIDs composing the OG main vessel that we're modifying (one of LAD, RCA, cfx)
			* intersection, 			dict {faceID: [ list of pointIDs in a branching vessel that also belong to main vessel]}
			* originals,				dict {faceID: np array of original point positions of inlets of that face}
			* affected_face_displace, 	dict {faceID: [list of displacement vectors (?) applied to this branch inlet]}
			* face_to_cap, 				dict {faceID: [list of pointIDs in cap at the end of the vessel labeled by faceID]}
			* point_connectivity, 		dict {pointID: [pointIDs that this point is in cell w/ ]}
			* easing, 					boolean, describes whether branch easing algorithm should be applied
			* boost, 					float, corresponding to scalar in front of branch shifting amount

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

		# get the mean displacement as the direction  
		average_displace = np.mean(displace_list, axis=0)
		# scale the mean displacement to preserve magnitude 
		average_displace /= np.linalg.norm(average_displace)
		average_displace *= np.linalg.norm(np.max(displace_list, axis=0))
		
		affected_face_displace[faceID] = average_displace
		
		#affected_face_displace[faceID] = np.mean(displace_list, axis=0)

	# apply displacement to affected branches and associated caps
	for faceID, displace in affected_face_displace.iteritems():
		cap_points = face_to_cap[faceID]
		vessel_points = face_to_points[faceID] - set(wall_region)
		branch_ids = cap_points.union(vessel_points)

		for pointID in branch_ids:
			cur_pt = wall_model.GetPoints().GetPoint(pointID)
			new_pt = [r + boost*dr for (r, dr) in zip(cur_pt, displace)]
			wall_model.GetPoints().SetPoint(pointID, new_pt)

		wall_model = branch_tilt(wall_model, intersection[faceID], originals[faceID], list(branch_ids))

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
	posts_zero_center = posts - np.mean(posts, axis=0)
	u_post,_,_ = np.linalg.svd(posts_zero_center.T)
	normal_post = u_post[-1]

	# determine how much to tilt 
	crossed = np.cross(normal_original, normal_post)
	G = np.array([
				[np.dot(normal_original, normal_post), -1*np.linalg.norm(crossed), 0],
				[np.linalg.norm(crossed), np.dot(normal_original, normal_post), 0 ], 
				[0, 0, 1]
				])

	vec_rejec = normal_post - (np.dot(normal_original, normal_post)*normal_original)
	vec_rejec /= np.linalg.norm(vec_rejec)

	Finv = np.zeros(G.shape)
	Finv[:,0] = normal_original
	Finv[:,1] = vec_rejec
	Finv[:,2] = -1*crossed

	tilt = np.matmul(Finv,G) 
	tilt = np.matmul(tilt, np.linalg.inv(Finv) )

	print '===================='
	print 'the tilt matrix is:'
	print tilt
	print 'normal original is:'
	print normal_original
	print 'normal post is:'
	print normal_post
	print 'testing tilt'
	print np.dot(tilt, normal_original)
	print 'F'
	print Finv
	print 'G'
	print G
	print '===================='

	# tilt branch 
	branch = extract_points(wall_model, pointIDs=branch_ids)
	branch -= np.mean(posts, axis=0)
	new_branch = np.dot(branch, tilt)
	new_branch += np.mean(posts, axis=0)

	# save tilts
	for pointID, point in zip(branch_ids, new_branch):
		wall_model.GetPoints().SetPoint(pointID, point)


	print 'done tilting branches'
	print '---------------------'

	return wall_model


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

	new_intersect = {face: [] for face in intersect.values()}
	for pt, face in intersect.iteritems(): new_intersect[face].append(pt)

	return (affected_face_displace, new_intersect)


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







