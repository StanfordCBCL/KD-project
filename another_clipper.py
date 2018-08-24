'''
	another_clipper.py
'''

import numpy as np
import vtk
from AneurysmGeneration.utils.slice import *



def clip_branch(polydata, normal, origin):
	'''
	'''

	plane = vtk.vtkPlane()
	plane.SetOrigin(origin)
	plane.SetNormal(-1*normal)

	extract = vtk.vtkExtractPolyDataGeometry()
	extract.SetInputData(polydata)
	extract.SetImplicitFunction(plane)
	extract.SetExtractBoundaryCells(True)
	extract.PassPointsOn()
	extract.Update()

	return extract.GetOutput()
	
def clip_branch_sphere(polydata, center, radius):
	'''
	'''

	sphere = vtk.vtkSphere()
	sphere.SetCenter(center)
	sphere.SetRadius(radius)

	extract = vtk.vtkExtractPolyDataGeometry()
	extract.SetInputData(polydata)
	extract.SetExtractInside(False)
	extract.SetImplicitFunction(sphere)
	extract.SetExtractBoundaryCells(True)
	extract.PassPointsOn()
	extract.Update()

	return extract.GetOutput()

def big_operator(suff, normals, origins, sphere_clips = None):

	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/LAD/ASI6/' + suff + '.vtp')
	reader.Update()
	poly_results = reader.GetOutput()

	
	for normal, origin in zip(normals, origins):
		poly_results = clip_branch(poly_results, normal, origin) 

	if sphere_clips is not None:
		print 'using sphere clips'
		for center, radius in zip(*sphere_clips):
			print center, radius
			poly_results = clip_branch_sphere(poly_results, center, radius)

	clipped_writer = vtk.vtkXMLPolyDataWriter()
	clipped_writer.SetInputData(poly_results)
	clipped_writer.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/LAD/' + suff + '.vtp')
	clipped_writer.Write()


if __name__ == "__main__":

	# data structure is going to be dict of {suff : list of lists containing normal, origin}

	suffs = ['lad1', 'lad2', 'lad3', 'lad4', 'lad5']

	# FOR LAD, ASI2
	# big_d = {
	# #'lad1': [np.array([[-12.14873847362136, 10.58989, -3.1367]]), np.array([[.63, .6467, -.4299]])], 
	# 'lad2': [np.array([[-12.1058, 10.605, -2.946], [-12.24216, 11.40157, -3.3217]]), np.array([[.4991, .7960, -.342366], [.63884, -.69974, .319761]])],
	# #'lad3': [np.array([[-12.04918, 10.51413, -2.9808], [-12.31338805, 11.38071, -3.53687]]), np.array([[.49559, .82418, -.27407], [-.337107, -.6744507, .65687]])],
	# #'lad4': [np.array([[-12.0736, 10.5614, -2.9405], [-12.29135, 11.42925, -3.60336]]), np.array([[.44635, .75312, -.483183], [-.2365302, -.676956, .696982]])],
	# #'lad5': [np.array([[-12.3028846, 10.6453, -2.9187], [-12.23315, 11.48406, -3.68404]]), np.array([[.433404, .851405, -.295415], [-.23567, -.69482, .6795]])]


	# }

	# FOR RCA, ASI6
	# big_d = {
	# 'p2': [np.array([[-6.123668, 13.9552, -7.0098]]), np.array([[-.300532, -.946231, .11882]])],
	# 'p3': [np.array([[-6.213015, 14.166392, -7.0079]]), np.array([[-.27811, -.95650, .08807]])],
	# 'p4': [np.array([[-6.2683, 14.29928, -6.73526]]), np.array([[-.275515, -.961023, .0229186]])],
	# 'p5': [np.array([[-5.93617, 14.29427, -6.8817]]), np.array([[-.262941, -.964801, .004606]])],
	# 'p5': [np.array([[-6.123668, 13.9552, -7.0098]]), np.array([[.300532, .946231, -.11882]])],


	# }

	# FOR LAD, ASI6
	big_d = {
	'lad1': [np.array([[-12.06955, 10.6508, -2.945514], [-13.7371712, 12.555614, -4.093051]]), np.array([[0.5519634, .775846, -.305639], [.8023378077, .400163, .442858423]])],

	# 'lad2': [np.array([[-11.98898, 10.57412, -3.0133259], [-13.99126, 12.571441, -3.99747] ]), np.array([[0.5644465, .65243, -.50571], [.945397, .22211, .23852]])],
	# 'lad3': [np.array([[-12.0323, 10.563, -3.1232], [-14.0746, 12.5727, -4.06032]]), np.array([[0.5467, .77858, -.30807], [.91761, .363155, .161603]])],
	# 'lad4': [np.array([[-12.0061, 10.626445, -2.958045], [-14.0991, 12.421, -4.029]]), np.array([[.563478, .756292, -.332438], [.91279, .37458, .162788]])], 
	# 'lad5': [np.array([[-11.9981, 10.601986, -2.9388431], [-14.1125, 12.4241255, -4.1772985]]), np.array([[.5359912, .81637837, -.2150345], [.91011496, .3682, .19006]])]
	

	}

	# sphere clips for lad asi6
	all_sphere_clips = {
	'lad1': [np.array([[-11.28277816, 12.2845688, -5.01866128]]), np.array([1.86445])],
	# 'lad2': [np.array([[-11.471234, 12.232, -4.846]]), np.array([-1.5641])],
	# 'lad3': [np.array([[-11.52133, 11.9257, -4.22575]]), np.array([.9483])],
	# 'lad4': [np.array([[-11.1808, 11.979, -4.76517]]), np.array([1.533202])],
	# 'lad5': [np.array([[-11.2788046, 11.9514549, -4.338]]), np.array([1.1758710])]

	}

	for suff, [origins, normals] in big_d.iteritems():
		print suff
		big_operator(suff, normals, origins, sphere_clips=all_sphere_clips[suff])
