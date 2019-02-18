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


fnames = ['p1', 'p5']
shapes = ['ASI4/', 'ASI6/']
read_dir = "clipped_results_short/RCA/"
out_path = "AdvectionDiffusion/"
total_path = "Artificial/RCA/"

for shape in shapes:
	for fname in fnames: 
		un_clipped = total_path + shape + fname + '/all_results.vtu'
		print 'the total number of nodes for ', un_clipped
		print return_unstructured(path).GetNumberOfPoints()

		gnids = extract_gnids(read_dir + shape + fname + '.vtu')
		with open(out_path + shape[:-1] + '_' + fname + '_nodes.txt', 'wb') as f: 
			gnids.tofile(f, sep="\n")
	
