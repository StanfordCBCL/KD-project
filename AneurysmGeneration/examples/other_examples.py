import numpy as np
import vtk


def gather_radii(wall_name, centerline=None):


	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall = wallreader.GetOutput()

	NoP = wall.GetNumberOfPoints()

	#preallocate np arrays

	radii = np.zeros(NoP)
	centroid = findCentroid(wall)

	print 'centroid is', centroid

	print 'NoP in wall is ', NoP

	for i in range(NoP):
		gNid = wall.GetPointData().GetArray('GlobalNodeID').GetTuple(i)
		pt = wall.GetPoints().GetPoint(i)
		d = vtk.vtkMath.Distance2BetweenPoints(pt, centroid)

		#print 'gNid is ', gNid
		#print 'pt is ', pt
		#print 'd is ', np.sqrt(d)

		radii[i] = np.sqrt(d)

	radii_vtk = nps.numpy_to_vtk(radii)
	radii_vtk.SetName('Radius')

	wall.GetPointData().AddArray(radii_vtk)

	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall)
	new.SetFileName('Results_' + wall_name)
	new.Write()
