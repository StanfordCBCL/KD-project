'''
visualization.py
'''

import vtk
from vtk.util import numpy_support as nps 
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_centerlines(centerlines, caplist, wall_name, cap_to_points)
	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

