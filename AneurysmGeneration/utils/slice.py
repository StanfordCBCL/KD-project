'''
	slice.py
'''

import numpy as np
import vtk
from parser import *
from pathreader import read_centerline
from normalization import normalized_centerline_pth
from batch import write_to_file

def determine_overlap(face_to_points, cap_to_points, NoP):
	'''

	input:

	output:


	'''

	print "determining overlap between walls"
	print "________________________________"

	point_to_face = {pointID:set() for pointID in range(NoP)}

	for faceID, pointIDs in face_to_points.iteritems():
		for faceID2, pointIDs2 in face_to_points.iteritems():
			intersect = pointIDs.intersection(pointIDs2) 
			if len(intersect) > 0:
				for pointID in intersect:
					point_to_face[pointID].add(faceID)
					point_to_face[pointID].add(faceID2)


	print "determining overlap between caps -> walls"
	print "________________________________"

	face_to_cap = {}
	for faceID, pointIDs in face_to_points.iteritems():
		for capID, pointIDs2 in cap_to_points.iteritems():
			intersect = pointIDs.intersection(pointIDs2)
			if len(intersect) > 0:
				face_to_cap[faceID] = pointIDs2
				break


	print "done determining overlap between walls and caps" 
	print "________________________________"

	return (point_to_face, face_to_cap)


def wall_isolation(face_list, cap_list, exclude, model_dir=None, wall_model=None, VALIDATION=True, EASING=False, PICKLE=False):
	'''

	input: 
		* face_list, a list of faces in the model (as ID #)
		* cap_list, a list of caps in the model (as ID #) 
		* exclude, a list of faces to exclude from the model (as ID #) 
		* model_dir, directory where the model files are
		* wall_name, name of the model wall 
		* VALIDATION, which should be included 

	output: 
		* dictionary that maps face_id to set of point IDs belonging to that face
	'''

	print 'computing structures from scratch: isolating wall sections'
	print '----------------------------------------------------------'

	# designate directories 
	if wall_model is None:
		wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"
		wallreader = vtk.vtkXMLPolyDataReader()
		wallreader.SetFileName(wall_name)
		wallreader.Update()
		wall_model = wallreader.GetOutput()

	if model_dir is None:
		model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"

	print "we will be preserving the following faces:"
	print face_list

	# get iteration boundaries
	NoP = wall_model.GetNumberOfPoints()
	NoC = wall_model.GetNumberOfCells()

	# initialize structures for holding preserved points corresponding to faces, caps
	face_to_points = {faceID:[] for faceID in face_list}
	cap_to_points = {capID:[] for capID in cap_list}


	point_connectivity = {}
	
	if EASING:
		print "the point connectivity data will be stored"


	for i in range(NoC):
		
		faceID = wall_model.GetCellData().GetArray('ModelFaceID').GetTuple(i)[0]
		cell_pt_ids = [int(wall_model.GetCell(i).GetPointIds().GetId(j)) for j in range(3)]

		if faceID in exclude:
			continue

		elif faceID in face_list:	
			face_to_points[faceID] += cell_pt_ids

			if EASING:
				# store the connectivity data
				for c, pointID in enumerate(cell_pt_ids):
					connected = [cell_pt_id for cell_pt_id in cell_pt_ids if cell_pt_id is not pointID]
					if pointID in point_connectivity.keys():
						point_connectivity[pointID] += connected
					else: 
						point_connectivity[pointID] = connected

		elif faceID in cap_list:
			cap_to_points[faceID] += cell_pt_ids
		

	# do a prelim processing 
	for faceID, pointIDs in face_to_points.iteritems():
		face_to_points[faceID] = set(pointIDs)

	for capID, pointIDs in cap_to_points.iteritems():
		cap_to_points[capID] = set(pointIDs)


	for pointID, connected in point_connectivity.iteritems():
		point_connectivity[pointID] = set(connected)

	if EASING:
		print 'done collecting point connectivity'
		print '----------------------------'




	print 'done isolating wall sections'
	print '----------------------------'

	if PICKLE:
		write_to_file("big_boy", (face_to_points, cap_to_points, point_connectivity, NoP))
	

	return (face_to_points, cap_to_points, point_connectivity, NoP)


def extract_points(polydata, pointIDs=None):
	'''
		Given an input polydata, extract points into ndarray of shape (NoP, 3)

		Optional pointIDs allows for return of only points corresponding to those pointIDs
	'''

	print 'extracting points'

	NoP = polydata.GetNumberOfPoints()
	points = np.zeros((NoP, 3))

	for i in range(NoP):
		points[i] = polydata.GetPoints().GetPoint(i)

	if pointIDs is None:
		return (NoP, points)
	else: 
		return points[pointIDs]

def obtain_expansion_region(wall_ref, NoP_wall, included_points, start=.1, end=.2, EPSILON=.01):
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

	wall_ref_axial = wall_ref[:, 0]
	wall_ref_theta = wall_ref[:, 1]

	# initialize datastructures
	wall_region_id = []
	start_id = []
	end_id = []

	# consider switching this for np.where in the future
	for i in included_points:
		if (wall_ref_axial[i] >= start) and (wall_ref_axial[i] <= end): 
			wall_region_id.append(i)

			if wall_ref_axial[i] < start + EPSILON:
				start_id.append(i) 

			if wall_ref_axial[i] > end - EPSILON:
				end_id.append(i)

	axial_pos = wall_ref[wall_region_id,0]
	theta_pos = wall_ref[wall_region_id,1]

	# report the number of border points identified:
	# print 'the number of points in the start border: ', len(start_border)
	# print 'the number of points in the end border: ', len(end_border)

	print 'Done obtaining the expansion region'
	print '-----------------------------------'

	return (wall_region_id, axial_pos, theta_pos, start_id, end_id) 


if __name__ == "__main__":

	print "testing slice.py"
	print "________________"

	wall_isolation()




