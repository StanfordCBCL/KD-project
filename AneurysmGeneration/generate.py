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
from utils.batch import *


def obtain_expansion_region(wall, centerline, included_points, start=.1, length=.1, EPSILON=.002):
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

	print 'Obtaining the expansion region'
	print '------------------------------'
	# project wall points to closest centerline point 
	normalized_wall, normalized_center, wall_to_center = projection(wall, centerline, included_points)

	# initialize datastructures
	wall_region_id = []
	axial_pos = []
	center_region_id = []	
	start_border = []
	end_border = []
	# determine how many points to iterate over
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center = len(centerline)

	# find the wall points projected into the desired centerline region
	for i in included_points:
		if (normalized_wall[i] >= start) and (normalized_wall[i] <= start + length): 
			wall_region_id.append(i)
			axial_pos.append(normalized_wall[i])
			if normalized_wall[i] < start + EPSILON:
				start_border.append(i)
			if normalized_wall[i] > start + length - EPSILON:
				end_border.append(i)

	# report the number of border poitns identified:
	print 'the number of points in the start border: ', len(start_border)
	print 'the number of points in the end border: ', len(end_border)

	# find the center points projected within the desired region
	for i in range(NoP_center):
		if (normalized_center[i] >= start) and (normalized_center[i] <= start + length):
			center_region_id.append(i)

	print 'Done obtaining the expansion region'
	print '-----------------------------------'
	return (wall_region_id, center_region_id, axial_pos, wall_to_center, start_border, end_border)


def shift_branches(wall_model, wall_region, intersection, affected_face_displace, face_to_cap, face_to_points, point_connectivity, easing):
	'''

	'''

	print 'Preparing to shift branches'
	print '------------------------------'

	# consolidate affected branch displacements
	for faceID, displace_list in affected_face_displace.iteritems():

		# compute the L2 norms of each displacement vector; select single displacement with the largest L2 norm
		displace_idx = np.argmax([np.sqrt(x**2 + y**2 + z**2) for (x,y,z) in displace_list])
		#average_displace = np.max(np.array(displace_list), axis=0)
		affected_face_displace[faceID] = displace_list[displace_idx]
		#affected_face_displace[faceID] = np.mean(displace_list, axis=0)

	# apply displacement to affected branches and associated caps
	for faceID, displace in affected_face_displace.iteritems():
		cap_points = face_to_cap[faceID]
		vessel_points = face_to_points[faceID] - set(wall_region)

		for pointID in cap_points.union(vessel_points):
			cur_pt = wall_model.GetPoints().GetPoint(pointID)
			new_pt = [r + dr for (r, dr) in zip(cur_pt, displace)]
			wall_model.GetPoints().SetPoint(pointID, new_pt)

		if easing:
			branch_easing(wall_model, intersection, vessel_points, point_connectivity)	

	print 'Done shifting branches'
	print '------------------------'
	return 

def branch_easing(wall_model, intersection, vessel_points, point_connectivity, num_iterations=5, aggress=.5):

	print 'Preparing to perform branch easing'
	print '----------------------------------'

	for i in range(num_iterations):

		print '>>>> Easing iteration # ', i

		easing = {}

		for pointID in intersection:
			connected = point_connectivity[pointID]

			for pointID_2 in connected:
				if pointID_2 in vessel_points:

					if pointID_2 in easing.keys():
						easing[pointID_2].append(pointID)
					else: 
						easing[pointID_2] = [pointID]


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


def acquire_start_end_radii(start_border, end_border, wall_model, wall_to_center):
	'''
		
	'''
	start_radii = acquire_radii(start_border, wall_model, wall_to_center)
	end_radii = acquire_radii(end_border, wall_model, wall_to_center)

	return (start_radii, end_radii)


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

def grow_aneurysm(wall_name, centerline, cur_face, face_to_points, point_to_face, face_to_cap, point_connectivity, options):
	'''
	input: 
		* name of wall vtp file
		* np array of centerline points 
		* included points, the pointIDs of the wall that we're considering
		* optional region specification
		* 

	output:
		* writes a file to current working directory named "modified_"+wall_name

	Given an input wall and centerline, artificially grow and aneurysm at the desired region. 
	'''

	print 'Preparing to grow aneurysm'
	print '--------------------------'


	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

	included_points = face_to_points[cur_face]

	(start, length, rad_max) = (options[k] for k in ['start', 'length', 'rad_max'])

	wall_region, center_region, axial_pos, wall_to_center, start_border, end_border = obtain_expansion_region(wall_model, 
																											centerline,
																											included_points,
																											start,
																											length)

	start_radii, end_radii = acquire_start_end_radii(start_border, end_border, wall_model, wall_to_center)
	
	affected_face_displace, intersect = organize_intersections(wall_region, point_to_face, cur_face)
	
	expand = interpolated_points(axial_pos, 
								(min(axial_pos), max(axial_pos)), 
								rad_shape=[np.median(start_radii), rad_max, np.median(end_radii) ] 
								)

	expand_np = np.zeros((1, wall_model.GetNumberOfPoints() ))


	for i, pointID in enumerate(wall_region):

		cur_pt = wall_model.GetPoints().GetPoint(pointID)
		normal = [r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[pointID]) ]

		displace = []
		new_pt = []

		if options['expansion_mode'] == 'scalar':
			displace = [expand[i]*dn for (r, dn) in zip(cur_pt, normal)]
			new_pt = [r + dr for (r, dr) in zip(cur_pt, displace)]
			

		elif options['expansion_mode'] == 'absolute':
			rad = np.linalg.norm(normal)
	 		normal_unit = [xi/rad for xi in normal]
			displace = [r*expand[i] for r in normal_unit]
			new_pt = [r + dr for (r, dr) in zip(wall_to_center[pointID], displace)]

			# after applying the displacement to the wall points, modify the magnitude of the displacement vector
			# to accomodate shifting of branches 
			displace = [d - n for (d, n) in zip (displace, normal)]


		# alter the current point's coordinates to reflect expansion
		wall_model.GetPoints().SetPoint(pointID, new_pt)

		# record the displacement for visualization on the np array
		expand_np[:,pointID] = np.linalg.norm(displace)

		# keep track of the displacement vectors affecting branch roots so that we know how much to shift
		# the branches
		if pointID in intersect.keys():
			affected_face_displace[intersect[pointID]].append(np.array(displace))

	# elif options['expansion_mode'] == 'absolute':

	# 	for i, pointID in enumerate(wall_region):

	# 		cur_pt = wall_model.GetPoints().GetPoint(pointID)
	# 		normal = [r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[pointID]) ]
	# 		rad = np.linalg.norm(normal)

	# 	 	normal_unit = [xi/rad for xi in normal]

	# 		displace = [r*expand[i] for r in normal_unit]
	# 		new_pt = [r + dr for (r, dr) in zip(wall_to_center[pointID], displace)]

	# 		wall_model.GetPoints().SetPoint(pointID, new_pt)

	# 		# record the displacement for visualization on the np array
	# 		expand_np[:,pointID] = np.linalg.norm(displace)

	# 		displace_adjusted = [d - n for (d, n) in zip (displace, normal)]
	# 		if pointID in intersect.keys():
	# 			affected_face_displace[intersect[pointID]].append(np.array(displace_adjusted))

	shift_branches(wall_model, wall_region, intersect.keys(), affected_face_displace, face_to_cap, face_to_points, point_connectivity, options['easing'])
	

	# add the expansion norms to the vtk file
	expand_vtk = nps.numpy_to_vtk(expand_np[0])
	expand_vtk.SetName('expand norm')
	wall_model.GetPointData().AddArray(expand_vtk)

	# write out the final vtk file
	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall_model)
	new.SetFileName(wall_name[:-4] + '_modified.vtp')
	new.Write()

	print 'done growing aneurysm!'
	print '----------------------'


def main():


	start = .4
	length = .3

	# define the location of models, centerlines, metadata and specify the wall_name
	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"

	EASING = False;

	options = aggregate_options(start=start, length=length, easing=EASING)
	print options

	# find the centerline files within the model directory and represent them as np arrays; 
	# find the names of the centerline files (without the .pth file ending)
	# note: this matches centerline name to the np array with all the point data
	centers, names = gather_centerlines(model_dir)

	# find the face IDs assigned to the cells in the model corresponding to the centerline names in the directory
	# note: this matches centerline name against its faceID 
	corresponding_faces, face_list = parse_facenames(names, model_dir)


	# be lazy -- hard code the exclude, face, and cap list
	exclude = [1]
	face_list = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
	cap_list = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

	# find the points corresponding to each relevant face ID and each cap ID
	# note: face_to_points is from faceID to list of pointID
	face_to_points, cap_to_points, point_connectivity, NoP = wall_isolation(face_list, cap_list, exclude, model_dir=model_dir, wall_name=wall_name, EASING=EASING)

	# determine the branching structure by looking at intersections of point IDs between face ID designations
	# determine which caps belong to which faces by looking at intersections
	point_to_face, face_to_cap = determine_overlap(face_to_points, cap_to_points, NoP)

	# prepare to grow an aneurysm at a specific region along a specified centerline
	# to do this, we input the set of centerline points as an np array of [xyz] and 
	# the set of pointIDs corresponding to the right wall region
	# 
	cur_name = names[2]
	cur_face = corresponding_faces[cur_name]
	cur_center = resample_centerline(centers[cur_name])
	cur_points = face_to_points[cur_face]
	grow_aneurysm(wall_name, cur_center, cur_face, face_to_points, point_to_face, face_to_cap, point_connectivity, options)


if __name__ == "__main__":

	main()

