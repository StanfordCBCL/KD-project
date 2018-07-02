'''
	generate.py

	Given an input wall vtp and centerline vtp, distort the wall to artificially generate
	an aneurysmal region. 

'''
import numpy as np
import vtk
from vtk.util import numpy_support as nps 

from mpl_toolkits.mplot3d import Axes3D
from utils.interpolation import *
from utils.normalization import *
from utils.parser import *
from utils.slice import *
from utils.batch import *
from utils.branch_ops import *


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
	wall_ref, normalized_center, wall_to_center, centerline_length = projection(wall, centerline, included_points)

	wall_ref_axial = wall_ref[:, 0]
	wall_ref_theta = wall_ref[:, 1]

	# compute the normalized length 
	length_frac = length/centerline_length

	# initialize datastructures
	wall_region_id = []
	axial_pos = []
	theta_pos = []
	center_region_id = []	
	start_border = []
	end_border = []
	# determine how many points to iterate over
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center = len(centerline)

	# find the wall points projected into the desired centerline region
	for i in included_points:
		if (wall_ref_axial[i] >= start) and (wall_ref_axial[i] <= start + length_frac): 
			wall_region_id.append(i)

			axial_pos.append(wall_ref_axial[i])
			theta_pos.append(wall_ref_theta[i])

			if wall_ref_axial[i] < start + EPSILON:
				start_border.append(i)
			if wall_ref_axial[i] > start + length_frac - EPSILON:
				end_border.append(i)

	# report the number of border points identified:
	print 'the number of points in the start border: ', len(start_border)
	print 'the number of points in the end border: ', len(end_border)

	# find the center points projected within the desired region
	for i in range(NoP_center):
		if (normalized_center[i] >= start) and (normalized_center[i] <= start + length_frac):
			center_region_id.append(i)

	print 'Done obtaining the expansion region'
	print '-----------------------------------'
	return (wall_region_id, center_region_id, axial_pos, theta_pos, wall_ref, wall_to_center, start_border, end_border)


def acquire_start_end_radii(start_border, end_border, wall_model, wall_to_center):
	'''
		
	'''
	start_radii = acquire_radii(start_border, wall_model, wall_to_center)
	end_radii = acquire_radii(end_border, wall_model, wall_to_center)

	return (start_radii, end_radii)



def grow_aneurysm(wall_name, face_to_points, point_to_face, face_to_cap, point_connectivity, cur_face, centerline, start, length, rad_max, easing, expansion_mode, suffix):

#def grow_aneurysm(wall_name, centerline, cur_face, face_to_points, point_to_face, face_to_cap, point_connectivity, options):
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

	wall_region, center_region, axial_pos, theta_pos, wall_ref, wall_to_center, start_border, end_border = obtain_expansion_region(
																											wall_model, 
																											centerline,
																											included_points,
																											start,
																											length)

	start_radii, end_radii = acquire_start_end_radii(start_border, end_border, wall_model, wall_to_center)
	
	affected_face_displace, intersect = organize_intersections(wall_region, point_to_face, cur_face)
	
	expand = interpolated_points(axial_pos, 
								(min(axial_pos), max(axial_pos)), 
								rad_shape=[np.mean(start_radii), rad_max, np.mean(end_radii) ] 
								)

	expand_np = np.zeros((1, wall_model.GetNumberOfPoints() ))


	for i, pointID in enumerate(wall_region):

		cur_pt = wall_model.GetPoints().GetPoint(pointID)
		wall_normal = [r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[pointID]) ]

		displace = []
		new_pt = []

		if expansion_mode == 'scalar':
			displace = [expand[i]*dn for (r, dn) in zip(cur_pt, normal)]
			new_pt = [r + dr for (r, dr) in zip(cur_pt, displace)]
			

		elif expansion_mode == 'absolute':

			# if pointID in start_border:
			# 	print 'skipping a start point'
			# 	continue

			rad = np.linalg.norm(wall_normal)
	 		wall_unit = [xi/rad for xi in wall_normal]
			displace = [r*expand[i] for r in wall_unit]
			new_pt = [r + dr for (r, dr) in zip(wall_to_center[pointID], displace)]

			# after applying the displacement to the wall points, modify the magnitude of the displacement vector
			# to accomodate shifting of branches 
			displace_adjusted = [d - n for (d, n) in zip (displace, wall_normal)]


		# alter the current point's coordinates to reflect expansion
		wall_model.GetPoints().SetPoint(pointID, new_pt)

		# record the displacement for visualization on the np array
		expand_np[:,pointID] = np.linalg.norm(displace)

		# keep track of the displacement vectors affecting branch roots so that we know how much to shift
		# the branches
		if pointID in intersect.keys():
			affected_face_displace[intersect[pointID]].append(np.array(displace_adjusted))


	print '>>>  the maximum displacement: ', np.max(expand_np)

	shift_branches(wall_model, wall_region, intersect.keys(), affected_face_displace, face_to_cap, face_to_points, point_connectivity, easing)
	

	# add the expansion norms to the vtk file
	expand_vtk = nps.numpy_to_vtk(expand_np[0])
	expand_vtk.SetName('expand norm')
	wall_model.GetPointData().AddArray(expand_vtk)

	wall_ref_transpose = np.transpose(wall_ref).copy()
	# add the normalized axial position wall array to the vtk file 
	norm_wall_vtk = nps.numpy_to_vtk(wall_ref_transpose[0, :])
	norm_wall_vtk.SetName('axial pos')
	wall_model.GetPointData().AddArray(norm_wall_vtk)

	# add the normalized theta position wall array to the vtk file
	theta_wall_vtk = nps.numpy_to_vtk(wall_ref_transpose[1, :])
	theta_wall_vtk.SetName('theta')
	wall_model.GetPointData().AddArray(theta_wall_vtk)

	# write out the final vtk file
	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall_model)
	new.SetFileName(wall_name[:-4] + '_modified_' + suffix + '.vtp')
	new.Write()


	print 'done growing aneurysm!'
	print '----------------------'


def main():


	# define the location of models, centerlines, metadata and specify the wall_name
	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"

	# some options 
	EASING = True
	PICKLE = False
	FROM_PICKLE = True
	BATCH = True
	PLOT_CL = True

	# find the centerline files within the model directory and represent them as np arrays; 
	# find the names of the centerline files (without the .pth file ending)
	# note: this matches centerline name to the np array with all the point data
	centers, names = gather_centerlines(model_dir)
	resampled = [resample_centerline(centers[name]) for name in names]


	# some random code for plotting the centerlines 
	'''
	if PLOT_CL:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		for centerline in resampled:
			ax.scatter(centerline[:,0], centerline[:,1], centerline[:,2])
		for name in names:
			center = centers[name]
			ax.scatter(center[:,0], center[:,1], center[:,2])
		plt.show()
	
	'''
	
	# find the face IDs assigned to the cells in the model corresponding to the centerline names in the directory
	# note: this matches centerline name against its faceID 
	corresponding_faces, face_list = parse_facenames(names, model_dir)


	# be lazy -- hard code the exclude, face, and cap list
	exclude = [1]
	face_list = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
	cap_list = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

	# find the points corresponding to each relevant face ID and each cap ID
	# note: face_to_points is from faceID to list of pointID
	face_to_points = None
	cap_to_points = None
	point_connectivity = None
	NoP = 0


	if FROM_PICKLE != True:
		print 'computing structures from scratch'
		print '---------------------------------'
		face_to_points, cap_to_points, point_connectivity, NoP = wall_isolation(face_list, cap_list, exclude, model_dir=model_dir, wall_name=wall_name, EASING=EASING)
	else:
		print 'reading structures from pickle'
		print '------------------------------'
		(face_to_points, cap_to_points, point_connectivity, NoP) = read_from_file("big_boy")

	if PICKLE:
		print 'writing structures to pickle'
		print '----------------------------'
		write_to_file("big_boy", (face_to_points, cap_to_points, point_connectivity, NoP))


	# determine the branching structure by looking at intersections of point IDs between face ID designations
	# determine which caps belong to which faces by looking at intersections
	point_to_face, face_to_cap = determine_overlap(face_to_points, cap_to_points, NoP)

	# prepare to grow an aneurysm at a specific region along a specified centerline
	# to do this, we input the set of centerline points as an np array of [xyz] and 
	# the set of pointIDs corresponding to the right wall region
	# 

	options = batch_targets(names, corresponding_faces, resampled, BATCH, EASING)
	print names
	print corresponding_faces

	for option in options: 

		print option

		grow_aneurysm(
			wall_name, 
			face_to_points, 
			point_to_face, 
			face_to_cap, 
			point_connectivity, 
			**option)

	
	
if __name__ == "__main__":

	main()

