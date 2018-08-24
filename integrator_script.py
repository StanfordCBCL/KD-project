'''
	integrator_script.py

	requires you to call this script by alias'ing pvpython in bashrc/bash_profile in order 
	to import paraview.simple module; :/ 


'''

#### import the simple module from the paraview
from paraview.simple import *
from vtk.util import numpy_support as nps	
import numpy as np

from AneurysmGeneration.utils.batch import write_to_file


def area_under_threshold(reader, tawss_upper=10.0, ):
	'''
		input:

		output:	

	'''

	threshold_bounds = np.arange(.05, tawss_upper, .5)
	area_fractions = np.zeros(threshold_bounds.shape)

	for i, curr_upper in enumerate(threshold_bounds):
		# create a new 'Threshold'
		threshold = Threshold(Input=reader)

		# Properties modified on threshold1
		threshold.Scalars = ['POINTS', 'vTAWSS']
		threshold.ThresholdRange = [0.0, curr_upper]

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



def main():

	computed_vars = {
	'WSS_THRESHOLD_AREA': {}, 
	'WSS_CYCLE': {},
	}

	corresponding_methods = {
	'WSS_THRESHOLD_AREA': area_under_threshold,
	'WSS_CYCLE': extract_wss_cycle
	}

	base_path = '/Users/alex/Documents/lab/KD-project/clipped_results_short/'
	# vessel = 'RCA/'
	vessel = 'LAD/'
	shapes = ['ASI2', 'ASI6']
	pos_sizes = None
	left_sizes = ['lad1', 'lad2', 'lad3', 'lad4', 'lad5']
	proximal = ['p1', 'p2', 'p3', 'p4', 'p5']
	medial = ['m1', 'm2', 'm3', 'm4', 'm5']
	distal = ['d1', 'd2', 'd3', 'd4', 'd5']

	for shape in shapes: 
		if 'RCA'in vessel: 
			pos_sizes = proximal + medial + distal
			if shape == 'ASI6': 
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
	