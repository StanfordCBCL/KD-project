from paraview.simple import *
import vtk

suff = 'LAD1'

# find source
vtp = FindSource('lad1.vtp')

# create a new 'Clip'
clip = Clip(Input=vtp)
clip.ClipType = 'Plane'

# init the 'Plane' selected for 'ClipType'
clip.ClipType.Origin = [-12.14873847362136, 10.58989, -3.1367]

# Properties modified on clip1
clip.Crinkleclip = 1

# Properties modified on clip1.ClipType
clip.ClipType.Normal = [.63, .6467, -.4299]

clip.UpdatePipeline()

clipped_writer = vtk.vtkXMLPolyDataWriter()
clipped_writer.SetInputData(clip.Output())
clipped_writer.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/' + suff + '.vtp')
#clipped_writer.Write()

