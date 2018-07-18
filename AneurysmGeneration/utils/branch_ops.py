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


def shift_branches(wall_model, wall_region, intersection, affected_face_displace, face_to_cap, face_to_points, point_connectivity, easing, boost=1.2):
	'''

	'''

	print 'Preparing to shift branches'
	print '------------------------------'

	# consolidate affected branch displacements
	for faceID, displace_list in affected_face_displace.iteritems():

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

		for pointID in cap_points.union(vessel_points):
			cur_pt = wall_model.GetPoints().GetPoint(pointID)
			new_pt = [r + boost*dr for (r, dr) in zip(cur_pt, displace)]
			wall_model.GetPoints().SetPoint(pointID, new_pt)

		if easing:
			branch_easing(wall_model, intersection, vessel_points, point_connectivity)	

	print 'Done shifting branches'
	print '------------------------'
	return 


def branch_easing(wall_model, intersection, vessel_points, point_connectivity, num_iterations=8, aggress=.8):

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

	intersect = {}
	affected_face_set = set()

	for pt in wall_region:
		belong_faces = point_to_face[pt]
		if len(belong_faces) > 1:
			for face in belong_faces: 
				if face != cur_face:
					intersect[pt] = face
					affected_face_set.add(face)

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







