"""baseline_extraction.py


"""

# external dependencies
import numpy as np
import vtk
from vtk.util import numpy_support as nps	
import sys
import os
import argparse

# internal modules
from AneurysmGeneration.utils.batch import *
from AneurysmGeneration.utils.normalization import *
from AneurysmGeneration.utils.slice import *



def main(): 
	args = parse_command_line(sys.argv)


if __name__ == "__main__": 
	main()