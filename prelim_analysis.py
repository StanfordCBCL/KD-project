'''
	prelim_analysis.py
'''

import numpy as np
import vtk
from vtk.util import numpy_support as nps
import sys
from AneurysmGeneration.utils.batch import *


def extract_points(polydata):
	'''
		Given an input polydata, extract points into ndarray of shape (NoP, 3)
	'''

	print 'extracting points'

	NoP = polydata.GetNumberOfPoints()
	points = np.zeros((NoP, 3))

	for i in range(NoP):
		points[i] = polydata.GetPoints().GetPoint(i)

	return points


def return_polydata(path):
	'''
		Given a file name, open and return the data
	'''
	print 'return polydata'

	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(path)
	reader.Update()
	polydata = reader.GetOutput()

	return polydata


def parse_command_line(args):
	
	print 'parsing command line'

	if len(args) == 1: 
		return ("AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp", "Artificial/RCA/p1/all_results.vtp")

	else: 
		return (args[1], args[2])

def apply_bounding_box(centerline, points, offsets=np.array([.07, .07, .07])):
	'''
		apply a bounding box based on min/max in axis 0 of centerline to points, 
		return the restricted indices


	'''

	print 'applying bb'

	upperbounds = np.amax(centerline, axis=0) + offsets
	lowerbounds = np.amin(centerline, axis=0) - offsets
	first = np.all(np.greater_equal(points, lowerbounds), axis=1)
	print first 

	inner = np.logical_and(
			np.all(np.greater_equal(points, lowerbounds), axis=1), 
			np.all(np.less_equal(points, upperbounds), axis=1) 
			)
	print inner
	patch_index = np.where(
		np.logical_and(
			np.all(np.greater_equal(points, lowerbounds), axis=1), 
			np.all(np.less_equal(points, upperbounds), axis=1) 
			) 
		)

	return patch_index


def min_dist(points_results, points_source):

	distances = np.sqrt(((points_results - points_source[:, np.newaxis])**2).sum(axis=2))
	return np.argmin(distances, axis=0)

def main():

	source_model_path, all_results_path = parse_command_line(sys.argv)

	# # we need the face to points correspondence from the original source model 
	# # we only need to get this once 
	(face_to_points, _, _, NoP) = read_from_file("big_boy")

	# # we need to get the source models
	# source_model_base_path = 'AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp'
	poly_source = return_polydata(source_model_path)

	# # we need the all_results_mesh
	# all_results_base_path = 'Artificial/RCA/p1/all_results.vtp'
	poly_results = return_polydata(all_results_path)
	
	gNid_array = nps.vtk_to_numpy(poly_results.GetPointData().GetArray('GlobalNodeID'))
	print 'the gnid array has shape', gNid_array.shape

	# # extract the raw points from the polydata
	points_source = extract_points(poly_source)
	points_results = extract_points(poly_results)

	# # load in the centerline
	centerline = read_from_file('RCA_cl')
	
	# use the centerline to apply a bounding box to the points that we need to look at 
	bounded_source_idx = apply_bounding_box(centerline, points_source)
	bounded_results_idx = apply_bounding_box(centerline, points_results)

	print points_results.shape
	print bounded_results_idx[0].shape

	print points_source.shape
	print bounded_source_idx[0].shape


	
	



	# # use np where to identify roughly close points? 

	# # 
	# # slim down the source points








if __name__ == "__main__":
	main()