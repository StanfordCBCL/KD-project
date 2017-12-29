import sys
import vtk
import numpy as np
import time
import os

# This function finds the centroid of the inputted vtk model
def findCentroid(input_model):

  numPts = input_model.GetNumberOfCells()

  centroid = [0.0, 0.0, 0.0]
  total_area = 0.0
  # Loop over all the surface triangles
  for icell in xrange(0, numPts):
  
    temp_cell = input_model.GetCell(icell)
    pts_cell = temp_cell.GetPointIds()
    
    # First, get the area of this cell
    vtkpt = temp_cell.GetPoints()
    p0 = vtkpt.GetPoint(0)
    p1 = vtkpt.GetPoint(1)
    p2 = vtkpt.GetPoint(2)
    temp_area = temp_cell.TriangleArea(p0,p1,p2)
    total_area = total_area + temp_area
    
    # Next, sum the coordinate values for each of the three nodal points
    x_sum = p0[0] + p1[0] + p2[0]
    y_sum = p0[1] + p1[1] + p2[1]
    z_sum = p0[2] + p1[2] + p2[2]
    
    # Now, compute the trapezoidal rule integration, multiply each summed
    # quantity by the area of the cell then divide by the number of vertices
    centroid[0] = centroid[0] + x_sum*temp_area/3.0
    centroid[1] = centroid[1] + y_sum*temp_area/3.0
    centroid[2] = centroid[2] + z_sum*temp_area/3.0
    
  centroid[0] = centroid[0]/total_area
  centroid[1] = centroid[1]/total_area
  centroid[2] = centroid[2]/total_area
  
  return centroid

if __name__ == '__main__':

  vtu_file = 'mesh-complete.mesh.vtu'
  
  # First, load in the vtu data from which we want to extract information
  vtu_reader = vtk.vtkXMLUnstructuredGridReader()
  vtu_reader.SetFileName(vtu_file)
  vtu_reader.Update()
  vtu_model = vtk.vtkUnstructuredGrid()
  vtu_model = vtu_reader.GetOutput()
  
  # Load in the inflow and outflow slices that we want to extract from the full model
  inflow_reader = vtk.vtkXMLPolyDataReader()
  inflow_reader.SetFileName(inflow_name)
  inflow_reader.Update()
  inflow_model = vtk.vtkPolyData()
  inflow_model = inflow_reader.GetOutput()
  inflow_numPts = inflow_model.GetNumberOfPoints()
  inflow_numCells = inflow_model.GetNumberOfCells()
  inflow_IDs = inflow_model.GetPointData().GetArray('GlobalNodeID')
  
  # Now, we want to "cut" a slice of the full model in the shape of the inlet and outlet slices.
  # To do this, I will first define a clip sphere that is centered on each of the inlets
  # and outlets with a radius that is big enough to encapsulate the entire inlet and outlet slices.
  # This cylinder cut will get rid of most of the full model. Then, I will define a slice
  # plane that is centered on the inlet and outlet slices with the same normal. This should
  # get the desired full model solution in the slice that I want.

  # Calculate the center of mass of the inlet slice by integrating each of the coordinates.
  # Since this mesh was generated using SimVascular, we assume that it is composed
  # completely of triangles. I am using trapezoidal rule here for integration
  centroid_inlet = findCentroid(inflow_model)
  print('Inlet centroid: (' + str(centroid_inlet[0]) + ', ' + str(centroid_inlet[1]) + ', ' + str(centroid_inlet[2]) + ')')
  
  # Now, calculate the area of the inlet slice
  masser = vtk.vtkMassProperties()
  masser.SetInputData(inflow_model)
  masser.Update()
  area_inlet = masser.GetSurfaceArea()
  
  # Use this area to calculate an approximate radius with which to perform a cylinder clip
  # We are going to inflate this a little bit since the inlet is not perfectly circular,
  # so we need a slightly larger area to make sure we do not clip around it
  radius_inlet = 1.5*np.sqrt(area_inlet/np.pi)
  print('Inlet radius sphere radius: ' + str(radius_inlet))

  # Now I want to perform a cut along the plane of the inlet slice. First need
  # the cell normals to the face
  normalGenerator = vtk.vtkPolyDataNormals()
  normalGenerator.SetInputData(inflow_model)
  normalGenerator.ComputePointNormalsOff()
  normalGenerator.ComputeCellNormalsOn()
  normalGenerator.Update()
  normals_test = normalGenerator.GetOutput()
  inlet_normal = normals_test.GetCellData().GetArray("Normals").GetTuple3(1)
  print("Inlet normal: (" + str(inlet_normal[0]) + ', ' + str(inlet_normal[1]) + ', ' + str(inlet_normal[2]) + ')')

  # Perform the cut
  cutPlane = vtk.vtkPlane()
  cutPlane.SetOrigin(centroid_inlet)
  cutPlane.SetNormal(inlet_normal)
  cutter = vtk.vtkCutter()
  cutter.SetCutFunction(cutPlane)
  cutter.SetInputData(vtu_model)
  cutter.Update()

  # Make the vtkSphere object
  inlet_sphere = vtk.vtkSphere()
  inlet_sphere.SetCenter(centroid_inlet)
  inlet_sphere.SetRadius(radius_inlet)

  inlet_clipper = vtk.vtkClipPolyData()
  inlet_clipper.SetClipFunction(inlet_sphere)
  inlet_clipper.InsideOutOn()
  inlet_clipper.SetInputData(cutter.GetOutput())
  inlet_clipper.Update()

  full_model_slice = inlet_clipper.GetOutput()
  
  # Write out the resulting inlet slice to verify
  inlet_slice_writer = vtk.vtkXMLPolyDataWriter()
  inlet_slice_writer.SetInputData(full_model_slice)
  inlet_slice_writer.SetFileName("inlet_slice_check.vtp")
  inlet_slice_writer.Write()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
