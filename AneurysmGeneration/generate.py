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


def grow_aneurysm(wall_name, face_to_points, point_to_face, face_to_cap, point_connectivity, cur_face, centerline, start, length, rad_max, easing, expansion_mode, suffix):
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

	
	# open the polydata
	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

	included_points = face_to_points[cur_face]
	NoP_wall, wall_points = extract_points(wall_model)

	#centerline=read_from_file('RCA_cl')

	wall_ref, normalized_center, wall_to_center, min_dists, centerline_length = projection(NoP_wall, centerline, wall_points, np.array(list(included_points)))

	# compute the normalized length -> normalized end position
	end=start+length/centerline_length

	wall_region, axial_pos, theta_pos, start_id, end_id = obtain_expansion_region(wall_ref, NoP_wall, included_points, start=start, end=end)

	start_border, end_border = wall_ref[start_id], wall_ref[end_id]
	start_radii, end_radii = min_dists[start_id], min_dists[end_id]
	# adjust = centerline_shift(start_id, end_id, axial_pos, wall_model, wall_to_center[start_id[0]], wall_to_center[end_id[-1]])
	# print adjust

	affected_face_displace, intersect = organize_intersections(wall_region, point_to_face, cur_face)

	originals = {faceID:extract_points(wall_model, pointIDs=intersect[faceID]) for faceID in intersect.keys()}
	
	# so we're actually going to clamp the end to the min radius so that the control flow condition later 
	# ensures some smoothness at the outlet 
	# turns out that this works really well! 
	expand = interpolated_points(axial_pos, 
								(min(axial_pos), max(axial_pos)), 
								rad_shape=[np.mean(start_radii), rad_max, np.min(end_radii) ] 
								)

	# expand = interpolation_2d(start_border, end_border, start_radii, end_radii, wall_ref[wall_region], start, end-start)

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

			if expand[i] < min_dists[pointID]: 
				expand[i] = min_dists[pointID]

	 		wall_unit = [xi/min_dists[pointID] for xi in wall_normal]
			displace = [r*expand[i] for r in wall_unit]
			new_pt = [r + dr for (r, dr) in zip(wall_to_center[pointID], displace)]

			# originally we try to use this centerline tilting thing but it's not that great, probably works a lot worse if we do it 
			# after applying the oriignla displacements; if we had used it initially to compute displacements, maybe that would be better. 
			# new_pt = [r - dr for (r, dr) in zip(new_pt, adjust[i])]

		# after applying the displacement to the wall points, modify displacement magnitude for branch shift
		displace_adjusted = [d - n for (d, n) in zip (displace, wall_normal)]

		for face, points in intersect.iteritems():
			if pointID in points:
# 				print '> recording disp for prop to this face:', face
				affected_face_displace[face].append(np.array(displace_adjusted))


		# alter the current point's coordinates to reflect expansion
		wall_model.GetPoints().SetPoint(pointID, new_pt)

		# record the displacement for visualization on the np array
		expand_np[:,pointID] = np.linalg.norm(displace)

	print '>>>  the maximum displacement: ', np.max(expand_np)

	shift_branches(wall_model, wall_region, intersect, originals, affected_face_displace, face_to_cap, face_to_points, point_connectivity, easing)
	

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
	EASING = False
	PICKLE = False
	FROM_PICKLE = True
	BATCH = True
	PLOT_CL = False

	# be lazy -- hard code the exclude, face, and cap list
	exclude = [1]
	face_list = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
	cap_list = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

	# find the centerline files within the model directory and represent them as np arrays; 
	# find the names of the centerline files (without the .pth file ending)
	# note: this matches centerline name to the np array with all the point data
	centers, names = gather_centerlines(model_dir)
	# resampled = [resample_centerline(centers[name]) for name in names]
	resampled = read_from_file('centerlines')
	
	# some random code for plotting the centerlines 
	if PLOT_CL: visualize_centerlines(names, centers, resampled)

	# find the face IDs assigned to the cells in the model corresponding to the centerline names in the directory
	# note: this matches centerline name against its faceID 
	corresponding_faces, face_list = parse_facenames(names, model_dir)

	print corresponding_faces

	# find the points corresponding to each relevant face ID and each cap ID
	# note: face_to_points is from faceID to list of pointID
	face_to_points = None
	cap_to_points = None
	point_connectivity = None
	NoP = 0

	if FROM_PICKLE != True:
		face_to_points, cap_to_points, point_connectivity, NoP = wall_isolation(face_list, cap_list, exclude, model_dir=model_dir, wall_model=wall_model, EASING=EASING, PICKLE=PICKLE)
	else:
		(face_to_points, cap_to_points, point_connectivity, NoP) = read_from_file("big_boy")


	# determine the branching structure by looking at intersections of point IDs between face ID designations
	# determine which caps belong to which faces by looking at intersections
	point_to_face, face_to_cap = determine_overlap(face_to_points, cap_to_points, NoP)

	options = batch_targets(names, corresponding_faces, resampled, BATCH, EASING)
	print names
	print corresponding_faces

	for option in options: 
		print option['suffix']
		grow_aneurysm(
			wall_name, 
			face_to_points, 
			point_to_face, 
			face_to_cap, 
			point_connectivity, 
			**option)

	
	
if __name__ == "__main__":

	main()

