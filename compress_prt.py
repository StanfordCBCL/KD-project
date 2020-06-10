"""compress_prt.py

A script for compressing advection-diffusion simulations so that there aren't 50000 files. 
Works by copying all of the fields of interest (e.g. "temperature") from many files into a single
combo .vtu file. 
"""
import vtk
import numpy as np
from vtk.util import numpy_support as nps
import argparse
import sys
from tqdm import tqdm

from AneurysmGeneration.utils.batch import return_unstructured 


def parse_commandline(argv): 
	"""Collects command line arguments passed in, returns them as a dictionary. Provides some usage documentation features. 
	
	Args:
	    argv (list): sys.argv for parsing using argparse
	
	Returns:
	    args (dict): dictionary mapping named commandline arguments to supplied values parsed out of sys.argv
	"""
	print 'parsing command line'

	parser = argparse.ArgumentParser(description='Usage for post processing the vtu result files of advection-diffusion simulations')
	parser.add_argument('--result_prefix', action="store", type=str, metavar='path_to/results_', required=True, help='path + prefix for files to be combined')
	parser.add_argument('--out_path', action="store", type=str, metavar='path_for/combo_file', required=True, help='output destination to write out the combined file')
	parser.add_argument('--start', action="store", type=int, required=True, help='first time step to be included in compressed file')
	parser.add_argument('--stop', action="store", type=int, required=True, help='last time step to be included in compressed file')
	parser.add_argument('--incr', action="store", type=int, default=50, help='number of time steps between output files to be compressed')
	
	args = vars(parser.parse_args())

	return args


def collect_field_of_interest(prefix, tsteps, extract_field='temperature', as_np=False):
	"""Collects the vtk point data array of interest from a series of separate .vtu files
	and returns a list of these arrays. Assumes that result files being parsed are named with format
	|prefix|+[time_step].vtu  
	
	Args:
	    prefix (str): path, naming prefix to results files that should be parsed
	    tsteps (list): list of integer time steps that should be included in the aggregated result file
	    extract_field (str, optional): name of the vtk point data array that should be extracted. by default, 
	        we will be extracting 'temperature' from advection-diffusion simulations.
	    as_np (bool, optional): optionally, convert vtk point data arrays to np arrays instead
	
	Returns:
	    compiled (list): list of collected arrays (1 from each of the result files), with names set to 
	        correspond to |extract_field| + [time_step]
	"""
	compiled = []
	for step in tqdm(tsteps, desc='extracting from results', file=sys.stdout): 
		results = return_unstructured(prefix + str(step) + '.vtu')
		field = results.GetPointData().GetArray(extract_field)
		field.SetName(extract_field+'_'+str(step))
		compiled.append(field)

	if not as_np: 
		return compiled
	else: 
		return [nps.vtk_to_numpy(field) for field in compiled]


def aggregate_fields_of_interest(first_file, output_dest, compiled, extract_field='temperature'): 
	"""Write the collected vtk point data arrays into a vtu file. 
	Arbitrarily picks one of the .vtu files to preserve points, connectivity, pressure, velocity
	information. 
	
	Args:
	    first_file (str): a str path to a vtu file that will serve as the template for copying in the 
	        remaining fields of interest.
	    output_dest (str): the filename for writing out a .vtu combined file
	    compiled (list): list of vtk point data arrays that should be copied into first_file
	    extract_field (str, optional): the name of the arrays that we're extracting 
	        this is used to remove the original array from the first file.
	
	"""
	dummy, reader = return_unstructured(first_file, return_reader=True)
	dummy.GetPointData().RemoveArray(extract_field)

	for field in tqdm(compiled, desc='writing results', file=sys.stdout): 
		dummy.GetPointData().AddArray(field)
		reader.Update()

	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetInputData(reader.GetOutput())
	writer.SetFileName(output_dest)
	writer.Write()


def main(): 

	args = parse_commandline(sys.argv)
	tsteps = list(range(args['start'], args['stop'], args['incr'])) + [args['stop']]

	compiled = collect_field_of_interest(args['result_prefix'], tsteps)
	aggregate_fields_of_interest(args['result_prefix'] + str(args['start']) + '.vtu', args['out_path'], compiled)


if __name__ == "__main__": 
	main()