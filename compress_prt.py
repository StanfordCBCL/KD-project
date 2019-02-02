"""compress_prt.py

A script for compressing advection-diffusion simulations so that there aren't 50000 files. 
"""
import vtk
import numpy as np
from vtk.util import numpy_support as nps
import argparse
import sys

from AneurysmGeneration.utils.batch import return_unstructured 


def parse_commandline(argv): 
	"""Summary
	
	Returns:
	    TYPE: Description
	"""
	print 'parsing command line'

	parser = argparse.ArgumentParser(description='Usage for post processing the vtu result files of advection-diffusion simulations')
	parser.add_argument('--result_prefix', action="store", type=str)
	parser.add_argument('--out', action="store", type=str)
	
	args = vars(parser.parse_args())

	return args


def collect_field_of_interest(prefix, extract_field='temperature', as_np=False):

	tsteps = ['000', '050'] + [str(num) + '0' for num in range(10, 105, 5)]
	compiled = []
	for step in tsteps: 
		results = return_unstructured(prefix + step + '.vtu')
		field = results.GetPointData().GetArray(extract_field)
		field.SetName(extract_field+'_'+step)
		compiled.append(field)

	if not as_np: 
		return compiled
	else: 
		return [nps.vtk_to_numpy(field) for field in compiled]


def combine(prefix, output_dest, compiled, extract_field='temperature'): 

	dummy, reader = return_unstructured(prefix + '000' + '.vtu', return_reader=True)
	dummy.GetPointData().RemoveArray(extract_field)

	for field in compiled: 
		dummy.GetPointData().AddArray(field)
		reader.Update()

	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetInputData(reader.GetOutput())
	writer.SetFileName(output_dest + '.vtu')
	writer.Write()


def main(): 

	args = parse_commandline(sys.argv)
	compiled = collect_field_of_interest(args['result_prefix'])
	combine(args['result_prefix'], args['out'], compiled)


if __name__ == "__main__": 
	main()