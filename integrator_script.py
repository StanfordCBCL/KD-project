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

def main():

	base_path = '/Users/alex/Documents/lab/KD-project/clipped_results_short/'
	vessel = 'RCA/'
	shapes = ['ASI2', 'ASI6']
	proximal = ['p1', 'p2', 'p3', 'p4', 'p5']
	medial = ['m1', 'm2', 'm3', 'm4', 'm5']
	distal = ['d1', 'd2', 'd3', 'd4', 'd5']

	all_area_fractions = {}

	for shape in shapes: 
		pos_sizes = proximal + medial + distal
		if shape == 'ASI6': pos_sizes = proximal
		for pos_size in pos_sizes: 
			full_path = base_path + vessel + shape + '/' + pos_size + '.vtp'
			print 'full_path is:', full_path

			reader = XMLPolyDataReader(FileName=[full_path])
			all_area_fractions[shape+'_'+pos_size] = area_under_threshold(reader)



	write_to_file('area_fractions', all_area_fractions)

if __name__ == "__main__":

	main()
	