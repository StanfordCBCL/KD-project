'''
	slice.py
'''

import numpy as np
import vtk
from parser import *
from pathreader import read_centerline
from normalization import normalized_centerline_pth



def wall_isolation(model_dir=None, wall_name=None, VALIDATION=True):
	'''

	input: 
		* directory where the model files are
		* name of the model wall 

	output: 
		* set of all points excluding the aorta
		* dictionary that maps face_id to set of point IDs belonging to that facee
	'''

	print "isolating wall sections"
	print "________________"

	# designate directories 
	if wall_name is None:
		wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"
	if model_dir is None:
		model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"

	# designate which face ids should be preserved
	good_faces = parse_facenames()
	face_to_points = {faceID:[] for faceID in good_faces}
	print good_faces

	# read in the wall
	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

	# get iteration boundaries
	NoP = wall_model.GetNumberOfPoints()
	NoC = wall_model.GetNumberOfCells()

	#record the IDs of points to be removed from consideration
	removed = []

	for i in range(NoC):
		
		faceID = wall_model.GetCellData().GetArray('ModelFaceID').GetTuple(i)[0]
		cell_pt_ids = wall_model.GetCell(i).GetPointIds()

		if faceID in good_faces:	
			for j in range(3):
				face_to_points[faceID].append(int(cell_pt_ids.GetId(j)))
		else:
			for j in range(3):
				removed.append(int(cell_pt_ids.GetId(j)))


	#perform set operations to remove duplicates 
	removed = set(removed)
	nonAortaPts = set(xrange(NoP)) - removed

	for faceID, pointIDs in face_to_points.iteritems():
		face_to_points[faceID] = set(pointIDs)

	if VALIDATION:
		# confirm that we haven't lost or created new points
		print len(removed), NoP, len(nonAortaPts)

		# confirm that our dictionary no longer contains duplicates
		c = 0
		for k, v in face_to_points.iteritems():
			print len(v)
			c += len(v)
		print c 

	return (nonAortaPts, face_to_points)

def face_grouping():





	'''
	mapped_center = np.zeros(NoP)
	centers = gather_centerlines(model_dir)

	for center in centers:	
		_, normalized, centerline_length = noramlized_centerline_pth(center)

	print centers.shape
	'''

if __name__ == "__main__":

	print "testing slice.py"
	print "________________"

	wall_isolation()




