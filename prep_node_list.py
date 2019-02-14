"""Summary
"""

# external dependencies
import numpy as np
import vtk
from vtk.util import numpy_support as nps	
import sys
import os
import argparse
from tqdm import tqdm

# internal modules
from AneurysmGeneration.utils.batch import *
from AneurysmGeneration.utils.normalization import *
from AneurysmGeneration.utils.slice import *


def extract_gnids(path): 
	unstructured_target = return_unstructured(path)
	return nps.vtk_to_numpy(unstructured_target.GetPointData().GetArray('GlobalNodeID'))


def get_num_nodes(path): 
	unstructured_target = return_unstructured(path)
	all_gnids = nps.vtk_to_numpy(unstructured_target.GetPointData().GetArray('GlobalNodeID'))
	return max(all_gnids)


def main(): 
	fnames = ['p1', 'p5']
	shapes = ['ASI2/', 'ASI4/', 'ASI6/']

	read_dir = "clipped_results_short/RCA/"
	out_path = "AdvectionDiffusion/"
	total_path = "Artificial/RCA/"
	max_nodes = []
	for shape in shapes:
		for fname in fnames: 
			gnids = extract_gnids(read_dir + shape + fname + '.vtu')
			with open(out_path + shape[:-1] + '_' + fname + '_nodes.txt', 'wb') as f: 
				for gnid in gnids: 
					f.write(str(gnid) + '\n')
		
			#max_nodes.append(get_num_nodes(total_path + shape + fname + '/all_results.vtu'))

	print max_nodes

if __name__ == "__main__": 
	main()