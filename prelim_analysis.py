'''
	prelim_analysis.py
'''

import numpy as np
import sys
from AneurysmGeneration.utils.batch import *
import argparse


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
	import vtk

	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(path)
	reader.Update()
	polydata = reader.GetOutput()

	return polydata


def parse_command_line(args):
	
	print 'parsing command line'

	parser = argparse.ArgumentParser(description='lol')
	parser.add_argument('--source', action="store", type=str, default="AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp")
	parser.add_argument('--results', action="store", type=str, default="Artificial/RCA/p1/all_results.vtp")
	parser.add_argument('--pkl', dest='feature', action='store_false')
	parser.add_argument('--vtk', dest='feature', action='store_true')
	parser.add_argument('--suff', action="store", type=str, default='p1')
	parser.set_defaults(feature=True)
	
	args = vars(parser.parse_args())

	return args

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


def prepare_chunks(idx, block_sz=100):
	'''
		input: 
			* idx, an np array of indices 
			* block_sz (optiona), the size of each partitioning block

	'''
	print 'preparing chunks'

	total = len(idx)
	print total
	partitions = block_sz*np.arange(total//block_sz)[1:]
	
	print partitions
	return partitions
	# return np.split(idx, partitions)


def min_dist(points_results, points_source):

	distances = np.sqrt(((points_results - points_source[:, np.newaxis])**2).sum(axis=2))
	return np.argmin(distances, axis=0)


def write_points_to_pkl(points_source, points_results, suff):
	write_to_file('points_' + suff, (points_source, points_results))


def main():


	CONTINUE=True
	args = parse_command_line(sys.argv)

	# # we need the face to points correspondence from the original source model 
	# # we only need to get this once 
	(face_to_points, _, _, NoP) = read_from_file("big_boy")

	points_source = None
	points_results = None

	if args['feature']: 
		source_model_path = args['source']
		all_results_path = args['results']

		# we need to get the source models
		# source_model_base_path = 'AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp'
		poly_source = return_polydata(source_model_path)

		# we need the all_results_mesh
		# all_results_base_path = 'Artificial/RCA/p1/all_results.vtp'
		poly_results = return_polydata(all_results_path)

		from vtk.util import numpy_support as nps	
		gNid_array = nps.vtk_to_numpy(poly_results.GetPointData().GetArray('GlobalNodeID'))

		print 'the gnid array has shape', gNid_array.shape

		# extract the raw points from the polydata
		points_source = extract_points(poly_source)
		points_results = extract_points(poly_results)

		write_points_to_pkl(points_source, points_results, args['suff'])

	else:
		'reading points from pickle'
		points_source, points_results = read_from_file('points_' + args['suff'])


	if CONTINUE:
		# # load in the centerline
		centerline = read_from_file('RCA_cl')
		
		# use the centerline to apply a bounding box to the points that we need to look at 
		bounded_source_idx = apply_bounding_box(centerline, points_source)[0]
		bounded_results_idx = apply_bounding_box(centerline, points_results)[0]

		print points_results.shape
		print bounded_results_idx.shape

		print points_source.shape
		print bounded_source_idx.shape

		# chunk one of the arrays so that we don't run out of memory when we 
		# perform the vectorized distance computation

		#partitions = prepare_chunks(bounded_results_idx) 
		mapping = np.zeros(points_results.shape[0])

		block_sz = 100
		for i in range(len(bounded_results_idx)//100):

			print 'split', i
			cur_idx = bounded_results_idx[block_sz*i:block_sz*(i+1)]
			mapping[cur_idx] = min_dist(points_results[cur_idx], points_source)

		write_to_file('mapping_'+args['suff'], mapping)

			

	
	



	# # use np where to identify roughly close points? 

	# # 
	# # slim down the source points








if __name__ == "__main__":
	main()