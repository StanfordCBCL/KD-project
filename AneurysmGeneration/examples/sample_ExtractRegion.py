    

meshfile = 'mymesh.vtp'

mesh= vtk.vtkXMLPolyDataReader()
mesh.SetFileName(meshfile)
mesh.GetOutput()
mesh.Update()

origin = [1,1,0]
normal = [1,0,0]

plane = vtk.vtkPlane()
plane.SetOrigin(origin[0] , origin[1], origin[2] )

plane.SetNormal(-normal[0], -normal[1],-normal[2])

clipper1=vtk.vtkExtractGeometry()
clipper1.SetImplicitFunction(plane)
clipper1.ExtractInsideOn()
clipper1.SetInputConnection(mesh.GetOutputPort())
clipper1.Update()
