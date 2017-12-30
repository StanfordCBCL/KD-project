import vtk 
from vtk.util import numpy_support as nps 
import numpy as np




def mutate_sphere(sphere_name, dr):

	sphere_vtk = vtk.vtkXMLPolyDataReader()
	sphere_vtk.SetFileName(sphere_name)
	sphere_vtk.Update()
	sphere = sphere_vtk.GetOutput()

	NoP = sphere.GetNumberOfPoints()

	new_pts = np.zeros((NoP, 3))

	centroid = findCentroid(sphere)

	pts = sphere.GetPoints()

	for i in range(NoP):

		prev_pt = sphere.GetPoints().GetPoint(i)
		newxyz = [r * (1 + dr) for r in prev_pt]
		new_pt = sphere.GetPoints().SetPoint(i, newxyz)

	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(sphere)
	new.SetFileName("modified_"+sphere_name)
	new.Write()


