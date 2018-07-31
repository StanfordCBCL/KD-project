'''
	prelim_analysis.py
'''

import numpy as np
import sys
import os
from AneurysmGeneration.utils.batch import *
from AneurysmGeneration.utils.normalization import *
from AneurysmGeneration.utils.slice import *
import argparse


def return_polydata(path, return_reader=False):
	'''
		Given a file name, open and return the data
	'''

	print 'return polydata'
	import vtk

	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName(path)
	reader.Update()
	polydata = reader.GetOutput()

	if return_reader:
		return (polydata, reader)

	return polydata


def parse_command_line(args):
	'''

	'''
	
	print 'parsing command line'

	parser = argparse.ArgumentParser(description='lol')
	parser.add_argument('--source', action="store", type=str, default="AneurysmGeneration/models/SKD0050/SKD0050_baseline_model_modified_p1.vtp")
	parser.add_argument('--results', action="store", type=str, default="Artificial/RCA/p1/all_results.vtp")
	parser.add_argument('--pkl', dest='pkl', action='store_true')
	parser.add_argument('--vtk', dest='vtk', action='store_true')
	parser.add_argument('--mapping', dest='mapping', action='store_true')
	parser.add_argument('--debug', dest='debug', action='store_true')
	parser.add_argument('--post', dest='post', action='store_true')
	parser.add_argument('--clip', dest='clip', action='store_true')

	parser.add_argument('--suff', action="store", type=str, default='p1')
	
	args = vars(parser.parse_args())

	return args


def apply_bounding_box(centerline, points, offsets=np.array([.11, .15, .8])):
	'''
		apply a bounding box based on min/max in axis 0 of centerline to points, 
		return the restricted indices
		Apply some offsets to each direction so that we definitely get all the points that we're looking for 

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

# def write_points_to_pkl(points_source, points_results, suff):
# 	write_to_file('points_' + suff, (points_source, points_results))

def main():

	args = parse_command_line(sys.argv)

	# we need the face to points correspondence from the original source model 
	# we only need to get this once 
	(face_to_points, _, _, NoP) = read_from_file("big_boy")

	points_source = None
	points_results = None

	if args['vtk']: 
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
		_, points_source = extract_points(poly_source)
		_, points_results = extract_points(poly_results)

		write_to_file('points_' + args['suff'], (points_source, points_results))

	elif args['pkl']:

		'reading points from pickle'
		points_source, points_results = read_from_file('points_' + args['suff'])


	if args['mapping']:

		# load in the centerline
		centerline = read_from_file('RCA_cl')
		
		# use the centerline to apply a bounding box to the points that we need to look at 
		bounded_source_idx = apply_bounding_box(centerline, points_source)[0]
		bounded_results_idx = apply_bounding_box(centerline, points_results)[0]

		# chunk one of the arrays so that we don't run out of memory when we 
		# perform the vectorized distance computation

		block_sz = 150
		n_splits = len(bounded_results_idx)//100

		mapping = np.zeros(points_results.shape[0])

		for i in range(n_splits):

			print '> split', i
			cur_idx = bounded_results_idx[block_sz*i:block_sz*(i+1)]

			if i == n_splits - 1:
				cur_idx = bounded_results_idx[block_sz*i:]

			mapping[cur_idx], _ = minimize_distances(points_results[cur_idx], points_source)

		write_to_file('mapping_'+args['suff'], mapping)

		vessel_ids = []
		mapped = np.zeros(mapping.shape[0])
		for c, cand in enumerate(mapping):
			if cand in face_to_points[8]:
				mapped[c] = 1
				vessel_ids.append(c)

		vessel_points = points_results[vessel_ids]

		write_to_file('mapped_'+args['suff'], (vessel_ids, vessel_points))

			
	if args['debug']:

		import vtk
		from vtk.util import numpy_support as nps

		mapping = read_from_file('mapping_'+args['suff'])

		print len(np.unique(mapping))
		print np.unique(np.isin(mapping, list(face_to_points[8])), return_counts=True)
		print np.where(
			np.isin(
				mapping, 
				face_to_points[8]
				)
			)
		mapped = np.zeros(mapping.shape[0])
		for c, cand in enumerate(mapping):
			if cand in face_to_points[8]:
				mapped[c] = 1
		# mapping[np.isin(mapping, face_to_points[8])] = 2
		mapping_vtk = nps.numpy_to_vtk(mapped)
		mapping_vtk.SetName('Mapped')

		all_results_path = args['results']
		poly_results = return_polydata(all_results_path)
		poly_results.GetPointData().AddArray(mapping_vtk)

		new=vtk.vtkXMLPolyDataWriter()
		new.SetInputData(poly_results)
		new.SetFileName('temp_'+ args['suff'] + '.vtp' )
		new.Write()

	
	if args['post']:
		import vtk
		from vtk.util import numpy_support as nps	
		all_results_path = args['results']
		poly_results, reader = return_polydata(all_results_path, return_reader=True)

		vessel_ids, vessel_points = read_from_file('mapped_'+args['suff'])
		# load in the centerline
		centerline = read_from_file('RCA_cl')

		# use normalization utils to map mesh points onto the centerline
		wall_ref, _, wall_to_center, min_dists, centerline_length = projection(NoP, centerline, points_results, vessel_ids)

		# get start and length from targets 
		targets = read_targets(as_dict=True)

		_, start, length, _ = targets[args['suff']]

		end=start+length/centerline_length

		wall_region, axial_pos, theta_pos, start_id, end_id = obtain_expansion_region(wall_ref, NoP, vessel_ids, start=start, end=end)
		
		v_tawss = nps.vtk_to_numpy(poly_results.GetPointData().GetArray('vTAWSS'))

		print np.mean(v_tawss[wall_region])

		if args['clip']:
			
			# get start and end point regions -- these are just start_id and end_id 

			points_start = extract_points(poly_results, start_id)
			points_end = extract_points(poly_results, end_id)

			# a single clip plane is defined by an origin and a normal
			# the origin can be computed as the average position of a ring of points 

			origin_start = np.mean(points_start, axis=0)
			origin_end = np.mean(points_end, axis=0)

			print origin_start

			# span, that we can use to ensure that the normals face in the right direction
			span = origin_end - origin_start
			span /= np.linalg.norm(span)

			# the normal can be computed as the average of the cross products of radial vectors
			# for vec in [points_start, points_end, origin_start]: 
			# 	vec.reshape(len(vec), 1)
			radial_start = points_start - origin_start
			radial_end = points_end - origin_end

			normal_start = np.cross(radial_start[0], radial_start[1])
			normal_end = np.cross(radial_end[0], radial_end[2])

			normal_start /= np.linalg.norm(normal_start)*np.sign(np.dot(normal_start, span))
			normal_end /= np.linalg.norm(normal_end)*np.sign(np.dot(normal_end, span))

			print normal_start
			# now, we have enough to build the clipper planes

			alpha = .02
			# shift the origin by a bit 
			origin_start -= span*alpha
			origin_end += span*alpha

			# define planes
			plane_start = vtk.vtkPlane()
			plane_start.SetOrigin(origin_start)
			plane_start.SetNormal(-1*span)

			plane_end = vtk.vtkPlane()
			plane_end.SetOrigin(origin_end)
			plane_end.SetNormal(span)

			# now let's pipe the two geometry extractors
			extract_start = vtk.vtkExtractPolyDataGeometry()
			extract_start.SetInputData(poly_results)
			extract_start.SetImplicitFunction(plane_start)
			extract_start.PassPointsOn()
			extract_start.Update()

			extract_end = vtk.vtkExtractPolyDataGeometry()
			extract_end.SetInputConnection(extract_start.GetOutputPort())
			extract_end.SetImplicitFunction(plane_end)
			extract_end.SetExtractBoundaryCells(True)
			extract_end.PassPointsOn()
			extract_end.Update()

			# we may have to; 
			# we can take origin = origin + (-ds)*n for a small ds in order to make sure that all the points from the origin
			# that we wanted to icnlud eare included 

			# next, we want to use a connectivity filter with SetExtractionModeToPointSeededRegions()
			connect = vtk.vtkPolyDataConnectivityFilter()
			connect.SetInputData(extract_end.GetOutput())
			connect.SetExtractionModeToPointSeededRegions()
			connect.AddSeed(wall_region[50]) # use an arbitrary point id from within the wall region to seed the connectivity filter
			connect.Update()

			region = connect.GetOutput()

			# now we write it out and pray that it worked
			clipped_writer = vtk.vtkXMLPolyDataWriter()
			clipped_writer.SetInputData(region)
			clipped_writer.SetFileName('clipped_results/RCA/'+args['suff']+'.vtp')
			clipped_writer.Write()



if __name__ == "__main__":
	main()



