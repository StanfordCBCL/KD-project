'''
	slice.py
'''

import numpy as np
import vtk
from parser import *
from pathreader import read_centerline
from normalization import noramlized_centerline_pth

def gather_centerlines(model_dir):
	centerlineNames = gather_centerline_names(model_dir)
	centers = []
	print model_dir
	print centerlineNames
	for centerline in centerlineNames:
		centers.append(read_centerline(centerline))

	return np.array(centers)


def slice_wall(model_dir):
	centers = gather_centerlines(model_dir)


	return None


if __name__ == "__main__":

	print "testing slice.py"
	print "________________"

	wall_name = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/SKD0050_baseline_model.vtp"
	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"

	good_faces = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

	NoP = wall_model.GetNumberOfPoints()
	NoC = wall_model.GetNumberOfCells()

	removed = []
	for i in range(NoC):
		
		faceID = wall_model.GetCellData().GetArray('ModelFaceID').GetTuple(i)[0]

		if faceID not in good_faces:
			cell_pt_ids = wall_model.GetCell(i).GetPointIds()
			for j in range(3):
				removed.append(int(cell_pt_ids.GetId(j)))


	removed = set(removed)
	nonAortaPts = set(xrange(NoP)) - removed
	print len(removed), NoP, len(nonAortaPts)

	mapped_center = np.zeros(NoP)
	centers = gather_centerlines(model_dir)

	for center in centers:
		_, normalized, centerline_length = noramlized_centerline_pth(center)

	print centers.shape
	#for i in nonAortaPts:
	#	pt = wall_model.GetPoints().GetPoint(i)




