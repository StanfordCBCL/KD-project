import sys
import vtk
import numpy as np
import time
import os

# ==============================================================================

# For each node in the mesh, finds the closest centerline point and records
# its normalized centerline coordinate
# 1. centerline_name - Name of the centerline file
# 2. mesh_model - Mesh of the vessel of interest
#      Will have a new data field called "Centerline_coordinate" which will
#      tell how far down the centerline each node in the mesh is
def find_normalized_coordinate(centerline_name, mesh_model):

  # Load in the centerline
  centerline_reader = vtk.vtkXMLPolyDataReader()
  centerline_reader.SetFileName(centerline_name)
  centerline_reader.Update()
  centerline_model = vtk.vtkPolyData()
  centerline_model = centerline_reader.GetOutput()
  centerline_numPts = centerline_model.GetNumberOfPoints()
  print("Number of points in the centerline: " + str(centerline_numPts))
  
  # Find the length of the centerline
  centerline_length = 0.0
  for i in xrange(1, centerline_numPts):
    pt = centerline_model.GetPoints().GetPoint(i)
    pt_prev = centerline_model.GetPoints().GetPoint(i-1)
    d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, pt_prev)
    d_temp = np.sqrt(d_temp)
    
    centerline_length = centerline_length + d_temp
  
  print("Centerline length: " + str(centerline_length))
  
  # Do another sweep and save the normalized coordinate for each point
  # on the centerline
  normalized_coordinates = vtk.vtkDoubleArray()
  normalized_coordinates.SetNumberOfComponents(1)
  normalized_coordinates.Allocate(centerline_numPts,10000)
  normalized_coordinates.SetNumberOfTuples(centerline_numPts)
  normalized_coordinates.SetName("Centerline_coordinate")
  normalized_coordinates.SetTuple1(0, 0.0)
  
  curr_distance = 0.0
  for i in xrange(1, centerline_numPts):
    pt = centerline_model.GetPoints().GetPoint(i)
    pt_prev = centerline_model.GetPoints().GetPoint(i-1)
    d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, pt_prev)
    d_temp = np.sqrt(d_temp)
    curr_distance = curr_distance + d_temp
    
    axial_coord = curr_distance / centerline_length
    normalized_coordinates.SetTuple1(i, axial_coord)
  
  centerline_model.GetPointData().AddArray(normalized_coordinates)
  writer = vtk.vtkXMLPolyDataWriter()
  writer.SetInputData(centerline_model)
  writer.SetFileName('centerline_with_coordinate.vtp')
  writer.Write()
  
  # Loop over each node in mesh_model and find the point on the centerline that
  # is closest to it. Assign a normalized coordinate according to this point
  model_numPts = mesh_model.GetNumberOfPoints()
  print("Number of points in mesh: " + str(model_numPts))
  projected_coordinate = vtk.vtkDoubleArray()
  projected_coordinate.SetNumberOfComponents(1)
  projected_coordinate.Allocate(model_numPts,10000)
  projected_coordinate.SetNumberOfTuples(model_numPts)
  projected_coordinate.SetName("Centerline_coordinate")
  
  for i in xrange(model_numPts):
    curr_node = mesh_model.GetPoint(i)
    min_distance = 999999.0
    min_index = -1
    
    # Loop over the centerline nodes to find the one that is closest
    for j in xrange(0, centerline_numPts):
      center_point = centerline_model.GetPoint(j)
      d_temp = vtk.vtkMath.Distance2BetweenPoints(curr_node, center_point)
      d_temp = np.sqrt(d_temp)
      if(d_temp < min_distance):
        min_distance = d_temp
        min_index = j
        
    coordinate_temp = normalized_coordinates.GetTuple1(min_index)
    projected_coordinate.SetTuple1(i, coordinate_temp)
    
  # Add the projected coordinate array to mesh_model
  mesh_model.GetPointData().AddArray(projected_coordinate)
  
  # Write out the model with centerline coordinate for verification
  datawriter = vtk.vtkXMLUnstructuredGridWriter()
  datawriter.SetInputData(mesh_model)
  
  
# ==============================================================================

if __name__ == '__main__':

  mesh_name = 'mesh-complete.mesh.vtu'
  center_name = 'lima_model_2_centerline.vtp'
  
  # Read in the .vtu mesh file
  datareader = vtk.vtkXMLUnstructuredGridReader()
  datareader.SetFileName(mesh_name)
  datareader.Update()
  mesh_model = vtk.vtkUnstructuredGrid()
  mesh_model = datareader.GetOutput()
  
  numNodes = mesh_model.GetNumberOfPoints()
  numCells = mesh_model.GetNumberOfCells()
  meshID = mesh_model.GetPointData().GetArray("GlobalNodeID")
  vtkNodes = mesh_model.GetPoints().GetData()
  vtkPointData = mesh_model.GetPointData()
  numpyNodes = vtk_to_numpy(vtkNodes)
  
  # Use find_normalized_coordinate to assign a normalized centerline coordinate
  # to each node in the mesh_model. After this function call, mesh_model should
  # have a new data entry called "Centerline_coordinate"
  find_normalized_coordinate(centerline_name, mesh_model)
