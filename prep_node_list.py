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


#fnames = ['p2', 'p3', 'p4']
fnames = ['p3']
shapes = ['ASI2/', 'ASI4/', 'ASI6/']
read_dir = "clipped_results_short/RCA/"
out_path = "AdvectionDiffusion/"
total_path = "Artificial/RCA/"

baseline_total = "Artificial/baseline/all_results.vtu"
baseline_base = "clipped_results_short/baseline/"
for shape in shapes:
	for fname in fnames: 
		# un_clipped = total_path + shape + fname + '/all_results.vtu'
		un_clipped = baseline_total
		print 'the total number of nodes for ', un_clipped
		print return_unstructured(un_clipped).GetNumberOfPoints()

		#clipped_name = read_dir + shape + fname + '.vtu'
		clipped_name = baseline_base + 'baseline_' + shape[:-1] + '_' + fname + '.vtu'
		gnids = extract_gnids(clipped_name)

		out_name = out_path + shape[:-1] + '_' + fname + '_nodes.txt'
		out_name = out_path + 'baseline_' + shape[:-1] + '_' + fname + '_nodes.txt'

		with open(out_name, 'wb') as f: 
			gnids.tofile(f, sep="\n")
	
