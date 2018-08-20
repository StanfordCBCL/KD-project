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
	


def big_operator(suff, normals, origins):

	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/RCA/ASI6/' + suff + '.vtp')
	reader.Update()
	poly_results = reader.GetOutput()

	
	for normal, origin in zip(normals, origins):
		poly_results = clip_branch(poly_results, normal, origin) 

	clipped_writer = vtk.vtkXMLPolyDataWriter()
	clipped_writer.SetInputData(poly_results)
	clipped_writer.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/RCA/' + suff + '.vtp')
	clipped_writer.Write()


if __name__ == "__main__":

	# data structure is going to be dict of {suff : list of lists containing normal, origin}

	suffs = ['lad1']#, 'lad2', 'lad3', 'lad4', 'lad5']

	# big_d = {
	# #'lad1': [np.array([[-12.14873847362136, 10.58989, -3.1367]]), np.array([[.63, .6467, -.4299]])], 
	# 'lad2': [np.array([[-12.1058, 10.605, -2.946], [-12.24216, 11.40157, -3.3217]]), np.array([[.4991, .7960, -.342366], [.63884, -.69974, .319761]])],
	# #'lad3': [np.array([[-12.04918, 10.51413, -2.9808], [-12.31338805, 11.38071, -3.53687]]), np.array([[.49559, .82418, -.27407], [-.337107, -.6744507, .65687]])],
	# #'lad4': [np.array([[-12.0736, 10.5614, -2.9405], [-12.29135, 11.42925, -3.60336]]), np.array([[.44635, .75312, -.483183], [-.2365302, -.676956, .696982]])],
	# #'lad5': [np.array([[-12.3028846, 10.6453, -2.9187], [-12.23315, 11.48406, -3.68404]]), np.array([[.433404, .851405, -.295415], [-.23567, -.69482, .6795]])]


	# }

	big_d = {
	'p2': [np.array([[-6.123668, 13.9552, -7.0098]]), np.array([[-.300532, -.946231, .11882]])],
	'p3': [np.array([[-6.213015, 14.166392, -7.0079]]), np.array([[-.27811, -.95650, .08807]])],
	'p4': [np.array([[-6.2683, 14.29928, -6.73526]]), np.array([[-.275515, -.961023, .0229186]])],
	'p5': [np.array([[-5.93617, 14.29427, -6.8817]]), np.array([[-.262941, -.964801, .004606]])],
	# 'p5': [np.array([[-6.123668, 13.9552, -7.0098]]), np.array([[.300532, .946231, -.11882]])],


	}

	for suff, [origins, normals] in big_d.iteritems():
		print suff
		big_operator(suff, normals, origins)
