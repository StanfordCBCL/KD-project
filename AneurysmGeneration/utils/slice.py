'''
	slice.py
'''

import numpy as np
import vtk
from parser import *
from pathreader import read_centerline
from normalization import normalized_centerline_pth


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


def wall_isolation(face_list, cap_list, exclude, model_dir=None, wall_name=None, VALIDATION=True, EASING=False):
	'''

	input: 
		* directory where the model files are
		* name of the model wall 

	output: 
		* dictionary that maps face_id to set of point IDs belonging to that face
	'''

	print "isolating wall sections"
	print "________________"

	# designate directories 
	if wall_name is None:
		wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"
	if model_dir is None:
		model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"

	print "we will be preserving the following faces:"
	print face_list

	# read in the wall
	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

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

	

	# if VALIDATION:
	# 	# confirm that we haven't lost or created new points
	# 	print len(removed), NoP, len(nonAortaPts)

	# 	# confirm that our dictionary no longer contains duplicates
	# 	c = 0
	# 	for k, v in face_to_points.iteritems():
	# 		print len(v)
	# 		c += len(v)
	# 	print c 

	print 'done isolating wall sections'
	print '----------------------------'

	return (face_to_points, cap_to_points, point_connectivity, NoP)



if __name__ == "__main__":

	print "testing slice.py"
	print "________________"

	wall_isolation()




