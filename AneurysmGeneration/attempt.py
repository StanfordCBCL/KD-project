import vtk
import numpy as np
from vtk.util import numpy_support as nps 
import os 
import matplotlib.pyplot as plt
from scipy import interpolate

from existing_scripts.clip_and_cut import *


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
	
def interpolated_points(x_interp, centerrange, rad_shape=None):

	if rad_shape is None:
		rad_shape = [0, .5, 1, .5, 0]

	print np.linspace(centerrange[0], centerrange[1], 3)
	print rad_shape

	tck = interpolate.splrep(np.linspace(centerrange[0], centerrange[1], 5), rad_shape, s=0)

	interpolated = interpolate.splev(x_interp, tck, der=0)
	# interpolated = np.interp(
	# 	x_interp, #np.linspace(centerrange[0], centerrange[1], numpoints), #(centerrange[1] - centerrange[0])/numpoints), 
	# 	np.linspace(centerrange[0], centerrange[1], 3),
	# 	rad_shape
	# 	)	

	return interpolated

def normalized_centerline(centerline_model):

	centerline_length = 0.0
	NoP = centerline_model.GetNumberOfPoints()
	normalized = {0.0: 0.0}

 	for i in range(1, NoP):
		pt = centerline_model.GetPoints().GetPoint(i)
		pt_prev = centerline_model.GetPoints().GetPoint(i-1)
		d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, pt_prev)
		d_temp = np.sqrt(d_temp)
		centerline_length = centerline_length + d_temp

	cur_length = 0.0
	for i in range(1, NoP):
		pt = centerline_model.GetPoints().GetPoint(i)
		pt_prev = centerline_model.GetPoints().GetPoint(i-1)
		d_temp = vtk.vtkMath.Distance2BetweenPoints(pt, pt_prev)
		d_temp = np.sqrt(d_temp)
		cur_length = cur_length + d_temp
		normalized[i] = cur_length/centerline_length


	return (NoP, normalized, centerline_length)

def projection(wall, centerline):
	NoP_wall = wall.GetNumberOfPoints()
	NoP_center, normalized_center, centerline_length = normalized_centerline(centerline)

	normalized_wall = np.zeros((NoP_wall))

	wall_to_center = {}

	for i in range(NoP_wall):
		wall_pt = wall.GetPoints().GetPoint(i)

		min_dist = float('inf')
		min_idx = -1
		for k in range(NoP_center):
			center_pt = centerline.GetPoints().GetPoint(k)
			cur_dist = vtk.vtkMath.Distance2BetweenPoints(wall_pt, center_pt)
			if cur_dist < min_dist:
				min_dist = cur_dist
				min_idx = k

		normalized_wall[i] = normalized_center[min_idx]

		wall_to_center[i] = centerline.GetPoints().GetPoint(min_idx)

	return (normalized_wall, normalized_center, wall_to_center)


def obtain_region(wall, centerline, start=.1, length=.1):
	normalized_wall, normalized_center, wall_to_center = projection(wall, centerline)

	wall_region_id = []
	axial_pos = []
	center_region_id = []	

	NoP_wall = wall.GetNumberOfPoints()
	NoP_center = centerline.GetNumberOfPoints()

	for i in range(NoP_wall):
		if (normalized_wall[i] >= start) and (normalized_wall[i] <= start + length): 
			wall_region_id.append(i)
			axial_pos.append(normalized_wall[i])


	for i in range(NoP_center):
		if (normalized_center[i] >= start) and (normalized_center[i] <= start + length):
			center_region_id.append(i)

	return (wall_region_id, center_region_id, axial_pos, wall_to_center)


def grow_aneurysm(wall_name, centerline_name, start=.1, length=.1):
	wallreader = vtk.vtkXMLPolyDataReader()
	wallreader.SetFileName(wall_name)
	wallreader.Update()
	wall_model = wallreader.GetOutput()

	centerreader = vtk.vtkXMLPolyDataReader()
	centerreader.SetFileName(centerline_name)
	centerreader.Update()
	centerline_model = centerreader.GetOutput()

	wall_region, center_region, axial_pos, wall_to_center = obtain_region(wall_model, centerline_model, start, length)

	#print axial_pos
	expand = interpolated_points(axial_pos, (min(axial_pos), max(axial_pos)) )

	#print len(expand)
	#print len(wall_region)
	#print len(axial_pos)

	for i, wall_id in enumerate(wall_region):
		#print 'setting point ', wall_id
		print 'expand is ', expand[i]
		cur_pt = wall_model.GetPoints().GetPoint(wall_id)
		normal = [r1 - r2 for (r1, r2) in zip(cur_pt, wall_to_center[wall_id]) ]
		print "cur", cur_pt
		print "n", normal
		new_pt = [r + expand[i]*dn for (r,dn) in zip(cur_pt, normal)]

		#[r * (1 + expand[i]) for r in cur_pt] #
		#print 'moving pt from ', cur_pt,
		#print 'to ', new_pt
		wall_model.GetPoints().SetPoint(wall_id, new_pt)

	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall_model)
	new.SetFileName("modified_"+wall_name)
	new.Write()



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



	return radii

	radii_vtk = nps.numpy_to_vtk(radii)
	radii_vtk.SetName('Radius')

	wall.GetPointData().AddArray(radii_vtk)

	new = vtk.vtkXMLPolyDataWriter()
	new.SetInputData(wall)
	new.SetFileName('Results_' + wall_name)
	new.Write()

	
def gen_cyl_and_center():
	

	points_list = [[0.0, 0.0, i] for i in np.linspace(-2.0, 5.0, 100)]


	# Create a vtkPoints object and store the points in it
	points = vtk.vtkPoints()

	for p in points_list:
		points.InsertNextPoint(p)

	centerpoly = vtk.vtkPolyData()

	centerpoly.SetPoints(points)

	writer2 = vtk.vtkXMLPolyDataWriter()
	writer2.SetFileName("test_model_centerline.vtp")
	writer2.SetInputData(centerpoly)
	writer2.Write()


def vis_radial_distribution(radii):
	return None



def main():

	#sphere_name = "Test.vtp"

	#mutate_sphere(sphere_name, 1)
	wall = "test_model.vtp"
	centerline = "test_model_centerline.vtp"

#	print interpolate(10, (0, 1))

	gen_cyl_and_center()

	grow_aneurysm(wall, centerline, start = .3, length=.1)

	#centerline = None
	#GatherRadii(wall)



if __name__ == "__main__":


	main()



