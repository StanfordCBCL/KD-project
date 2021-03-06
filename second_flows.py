"""
	second_flows.py

	on my Mac, I was able to only able to import paraview.simple by by alias'ing pvpython in bashrc/bash_profile 
	otherwise, this should work pretty directly


"""

# external dependencies
from paraview.simple import * # import the simple module from the paraview
import numpy as np
import vtk
from vtk.util import numpy_support as nps	


# internal dependencies
from AneurysmGeneration.utils.batch import *
# from AneurysmGeneration.utils.normalization import *
# from AneurysmGeneration.utils.slice import *
from AneurysmGeneration.utils.parser import *
# from AneurysmGeneration.utils.slice import determine_overlap

# get the normal, centerline to define the slice  <<<<
# do the slice an dget the output <<<<
# take the output and do the integrand computation
# write a new array???  <<<< 
# do integrate variables 


def compute_secondary_flow(centers, normals, suffix): 

	filename = '/Users/alex/Documents/lab/KD-project/clipped_results_short/RCA/ASI2/' + suffix + '.vtu'
	vtu = XMLUnstructuredGridReader(FileName=[filename])

	step_lower = 3000
	step_upper = 4000
	tstep = 50
	steps = np.arange(step_lower, step_upper + tstep, tstep)

	all_second_flows = []
	areas = []
	velocity_peaks = []

	for center, normal in zip(centers, normals):

		second_flows_slice = np.zeros(len(steps))

		cross_slice = Slice(Input=vtu)
		cross_slice.SliceType = 'Plane'
		cross_slice.SliceOffsetValues = [0.0]

		for i, step in enumerate(steps):

			cross_slice.SliceType.Origin = center
			cross_slice.SliceType.Normal = normal

			## use the sliced field result to manually subtract the out-of-plane component 
			## this is when we compute u - (u \cdot n)n 
			slice_results = paraview.servermanager.Fetch(cross_slice)
			sliced_field = nps.vtk_to_numpy(slice_results.GetPointData().GetArray('velocity_0' + str(step)))
			
			inplane_manual = sliced_field - normal * np.dot(sliced_field, normal)[:, None]

			## (u - (u \cdot n)n )^2  -- chose to compute this as squared L2 norm because u - (u \cdot n)n is a vector
			integrand = np.square(np.linalg.norm(inplane_manual, axis=1)) 

			integrand_vtk = nps.numpy_to_vtk(integrand)
			integrand_vtk.SetName('integrand_0' + str(step))

			slice_results.GetPointData().AddArray(integrand_vtk)
			slice_results.GetPointData().Update()

			temp = vtk.vtkXMLPolyDataWriter()
			temp.SetInputData(slice_results)
			temp.SetFileName('cache.vtu' )
			temp.Write()

			cache_slice = XMLPolyDataReader(FileName=['cache.vtu'])

			integrateVariables = IntegrateVariables(Input=cache_slice)
			integrated_result = paraview.servermanager.Fetch(integrateVariables)

			second_flows_slice[i] = integrated_result.GetPointData().GetArray('integrand_0' + str(step)).GetTuple(0)[0]
		
		all_second_flows.append(second_flows_slice)


		integrateVariables = IntegrateVariables(Input=cache_slice)
		integrated_result = paraview.servermanager.Fetch(integrateVariables)

		areas.append(integrated_result.GetCellData().GetArray('Area').GetTuple(0)[0])

		integrated_velocity_peak = max(map(np.linalg.norm, [
			integrated_result.GetPointData().GetArray('velocity_03200').GetTuple(0),
			integrated_result.GetPointData().GetArray('velocity_03250').GetTuple(0), 
			integrated_result.GetPointData().GetArray('velocity_03300').GetTuple(0), 
			integrated_result.GetPointData().GetArray('velocity_03350').GetTuple(0), 

		]))

		velocity_peaks.append(integrated_velocity_peak)

	return all_second_flows, areas, velocity_peaks


def normalized_centerline_pth(center):
	'''
	input:
		* np array of shape (NoP, 3)

	output: 
		* NoP
		* np array of length NoP, containing normalized coordinate for each point
		* total centerline length

	Assigns each centerline point a total length-normalized position, holding assigned coordinate
	in form of np array with shape (NoP,). 
	
	'''

	print 'normalizing the centerline'
	print '--------------------------'
	centerline_length = 0.0
	NoP = len(center)
	normalized = np.zeros(NoP)

	for i in range(1, NoP):
		pt = center[i]
		prev_pt = center[i-1]
		d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, prev_pt)
		d_temp = np.sqrt(d_temp)
		centerline_length += d_temp
		normalized[i] = centerline_length

	normalized /= centerline_length

	return (NoP, normalized, centerline_length)


def slice_parameters(centerline, start, length, positions=6): 

	## watch!! 
	centerline=read_from_file('RCA_cl')

	NoP_center, normalized_center, centerline_length = normalized_centerline_pth(centerline)

	# compute the normalized length -> normalized end position
	end = start + length/centerline_length

	# print(len(centerline))
	aneurysm_segment = centerline[(normalized_center >= start) & (normalized_center <= end)]

	NoP = aneurysm_segment.shape[0]
	p0 = np.roll(aneurysm_segment, shift = -1, axis= 0)
	p1 = np.roll(aneurysm_segment, shift = 0, axis = 0)
	p2 = np.roll(aneurysm_segment, shift = 1, axis = 0)

	t21 = p2 - p1
	t21_normed = np.divide(t21, np.linalg.norm(t21, axis=1).reshape(NoP, 1))

	# print(len(aneurysm_segment))
	
	interval = len(aneurysm_segment)/(positions+1)
	indices = [i*interval for i in range(1, positions+1)]

	centers = aneurysm_segment[indices]
	normals = t21_normed[indices]

	return centers, normals


def read_batch_parameters(names, resampled, targets_name):
	'''

	'''

	targets = read_targets(targets_name)
	
	options = []
	for target in targets:
		(vessel, start, length, rad_max, suffix) = target
		options.append({
			'start': start,
			'length': length,
			'rad_max': rad_max,
			'centerline': resampled[vessel], 
			'suffix': suffix})

	return options


def batch_second_flows(): 

	# define the location of models, centerlines, metadata and specify the wall_name
	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	targets_name = "/AneurysmGeneration/targets.txt"
	# targets_name = "/AneurysmGeneration/left_targets.txt"

	# find the centerline files within the model directory and represent them as np arrays; 
	# find the names of the centerline files (without the .pth file ending)
	# note: this matches centerline name to the np array with all the point data
	_, names = gather_centerlines(model_dir)
	resampled = read_from_file('centerlines')

	options = read_batch_parameters(names, resampled, targets_name)	

	results = {}
	for option in options: 
		print option['suffix']
		centers, normals = slice_parameters(option['centerline'], option['start'], option['length'])
		all_second_flows, areas, velocity_peaks = compute_secondary_flow(centers, normals, option['suffix'])

		results[option['suffix']] = (all_second_flows, areas, velocity_peaks, centers, normals)

	write_to_file('second_flows_dict_2', results)


if __name__ == "__main__": 
	# main()
	batch_second_flows()


