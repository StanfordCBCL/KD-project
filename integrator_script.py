'''
	integrator_script.py

	on my Mac, I was able to only able to import paraview.simple by by alias'ing pvpython in bashrc/bash_profile 
	otherwise, this should work pretty directly

	script function: 
		define data structures for holding integrated results 

		define some strings used to navigate the data directory 
		assumes that the data files are located in directories structured like: 
			[base_path] 
				[vessel] 
					[shape_index] 
						[position, size]

		the integrated variable dictionaries are written out as pickles


'''

# external dependencies
from paraview.simple import * # import the simple module from the paraview
from vtk.util import numpy_support as nps	
import numpy as np

# internal dependencies
from AneurysmGeneration.utils.batch import write_to_file


def area_under_threshold(reader, tawss_upper=10.0):
	'''
	
	Args:
	    reader (XMLPolyDataReader): XMLPolyDataReader object used by pv simple functions to threshold, integrate contents of polydata file. 
	    tawss_upper (float, optional): optional upper bound for evaluating the time-averaged wall shear stress. 
	
	Returns:
	    area_fractions (nd_array): fractional surface areas exposed to successively increasing thresholds of tawss
	
	'''

	threshold_bounds = np.arange(.0, tawss_upper, .5)
	threshold_bounds[0] += .05
	area_fractions = np.zeros(threshold_bounds.shape)

	for i, curr_upper in enumerate(threshold_bounds):
		threshold = Threshold(Input=reader)

		threshold.Scalars = ['POINTS', 'vTAWSS']
		threshold.ThresholdRange = [0.0, curr_upper]

		integrateVariables = IntegrateVariables(Input=threshold)
		result = paraview.servermanager.Fetch(integrateVariables)

		try:
			area_fractions[i] = result.GetCellData().GetArray('Area').GetTuple(0)[0] 
		except:
			area_fractions[i] = 0

	integrateVariables = IntegrateVariables(Input=reader)
	result = paraview.servermanager.Fetch(integrateVariables)
	total_area = result.GetCellData().GetArray('Area').GetTuple(0)[0]

	area_fractions /= total_area

	print 'the area fractions under threshold for TAWSS w/ upper bound: ', tawss_upper
	print area_fractions

	return area_fractions


def extract_wss_cycle(reader, step_lower=3000, step_upper = 4000, tstep = 50):
	'''

		compute the average wall shear stress at every point in the cardiac cycle 
		compute the average TAWSS 
	'''

	steps = np.arange(step_lower, step_upper + tstep, tstep)
	wss_cycle_vec = np.zeros((len(steps), 3))

	for i, step in enumerate(steps):
		integrateVariables = IntegrateVariables(Input=reader)
		result = paraview.servermanager.Fetch(integrateVariables)

		wss_cycle_vec[i] = result.GetPointData().GetArray('vWSS_0' + str(step)).GetTuple(0)

	total_area = result.GetCellData().GetArray('Area').GetTuple(0)[0]

	wss_cycle = np.linalg.norm(wss_cycle_vec, axis=1)/total_area
	avg_vtawss = np.linalg.norm(result.GetPointData().GetArray('vTAWSS').GetTuple(0))/total_area

	return (wss_cycle, avg_vtawss)


def osi_above_threshold(reader, osi_upper=.5):
	'''
		input:

		output:	

	'''

	threshold_bounds = np.linspace(0, osi_upper, num=75)
	area_fractions = np.zeros(threshold_bounds.shape)

	for i, cur_lower in enumerate(threshold_bounds):
		# create a new 'Threshold'
		threshold = Threshold(Input=reader)

		# Properties modified on threshold1
		threshold.Scalars = ['POINTS', 'vOSI']
		threshold.ThresholdRange = [cur_lower, 10.0]

		# create a new 'Integrate Variables'
		integrateVariables = IntegrateVariables(Input=threshold)
		result = paraview.servermanager.Fetch(integrateVariables)

		try:
			area_fractions[i] = result.GetCellData().GetArray('Area').GetTuple(0)[0] #unwrap this shit 
		except:
			area_fractions[i] = 0

	# compute the total area
	integrateVariables = IntegrateVariables(Input=reader)
	result = paraview.servermanager.Fetch(integrateVariables)
	total_area = result.GetCellData().GetArray('Area').GetTuple(0)[0]

	area_fractions /= total_area
	print area_fractions

	return area_fractions


def main():
	"""
	"""

	# define some data structures in which the integrated results are going to be held
	computed_vars = {
	'WSS_THRESHOLD_AREA': {}, 
	'WSS_CYCLE': {},
	'OSI_THRESHOLD_AREA':{}
	}

	corresponding_methods = {
	'WSS_THRESHOLD_AREA': area_under_threshold,
	'WSS_CYCLE': extract_wss_cycle,
	'OSI_THRESHOLD_AREA': osi_above_threshold
	}


	# define some parameters for finding the data .vtp files 
	base_path = '/Users/alex/Documents/lab/KD-project/clipped_results_short/'
	vessel = 'RCA/'
#	vessel = 'LAD/'
#	shapes = ['ASI2', 'ASI6']
	shapes = ['ASI2', 'ASI4', 'ASI6']
	pos_sizes = None
	left_sizes = ['lad1', 'lad2', 'lad3', 'lad4', 'lad5']
	proximal = ['p1', 'p2', 'p3', 'p4', 'p5']
	medial = ['m1', 'm2', 'm3', 'm4', 'm5']
	distal = ['d1', 'd2', 'd3', 'd4', 'd5']

	# iterate over the data directory 
	for shape in shapes: 
		if 'RCA'in vessel: 
			pos_sizes = proximal + medial + distal
			if shape == 'ASI6' or shape == 'ASI4': 
				pos_sizes = proximal
				
		if 'LAD' in vessel: 
			pos_sizes = left_sizes

		for pos_size in pos_sizes: 
			full_path = base_path + vessel + shape + '/' + pos_size + '.vtp'
			print 'working on:', full_path

			reader = XMLPolyDataReader(FileName=[full_path])

			for var, val in computed_vars.iteritems():
				if val is not None:
					computed_vars[var][shape+'_'+pos_size] = corresponding_methods[var](reader)


	for var, val in computed_vars.iteritems():
		if val is not None: write_to_file(var+'_' + vessel[:-1], val)

if __name__ == "__main__":
	main()
	