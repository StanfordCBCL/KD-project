'''
	another_clipper.py
'''

import numpy as np
import vtk
from AneurysmGeneration.utils.slice import *
from AneurysmGeneration.utils.batch import return_unstructured, return_polydata



def clip_branch(polydata, normal, origin, extractor):
	'''
	'''

	plane = vtk.vtkPlane()
	plane.SetOrigin(origin)
	plane.SetNormal(-1*normal)

	extract = extractor
	extract.SetInputData(polydata)
	extract.SetImplicitFunction(plane)
	extract.SetExtractBoundaryCells(True)
#	extract.PassPointsOn()
	extract.Update()

	return extract.GetOutput()

	
def clip_branch_sphere(polydata, center, radius, extractor):
	'''
	'''

	sphere = vtk.vtkSphere()
	sphere.SetCenter(center)
	sphere.SetRadius(radius)

	extract = extractor
	extract.SetInputData(polydata)
	extract.SetExtractInside(False)
	extract.SetImplicitFunction(sphere)
	extract.SetExtractBoundaryCells(True)
#	extract.PassPointsOn()
	extract.Update()

	return extract.GetOutput()


def apply_clipping(path, vessel, shape, suff, 
					normals, origins, 
					reader, writer, extractor, extension,
					sphere_clips = None): 
	"""Summary
	
	Args:
	    path (TYPE): Description
	    vessel (TYPE): Description
	    shape (TYPE): Description
	    suff (TYPE): Description
	    normals (TYPE): Description
	    origins (TYPE): Description
	    reader (TYPE): Description
	    writer (TYPE): Description
	    extractor (TYPE): Description
	    extension (TYPE): Description
	    sphere_clips (None, optional): Description
	"""
	print 'entering apply_clipping'

	results = reader(path + vessel + shape + suff + extension)

	for normal, origin in zip(normals, origins): 
		results = clip_branch(results, normal, origin, extractor)

	if sphere_clips is not None: 
		print 'using sphere clips'
		for center, radius in zip(*sphere_clips): 
			print center, radius
			results = clip_branch_sphere(results, center, radius, extractor)

	clip_writer = writer
	clip_writer.SetInputData(results)
	clip_writer.SetFileName(path + vessel + shape + suff + extension)
	clip_writer.Write()

	print 'exiting apply_clipping'


def big_operator(suff, normals, origins, sphere_clips = None):

	reader = vtk.vtkXMLPolyDataReader()
	reader.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/RCA/ASI4/' + suff + '.vtp')
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
	clipped_writer.SetFileName('/Users/alex/Documents/lab/KD-project/clipped_results_short/RCA/' + suff + '.vtp')
	clipped_writer.Write()


def get_clip_parameters(vessel, shape): 

	# FOR RCA, ASI4
	rca_asi4 = {
	'p5': [np.array([[-6.182215893137097, 14.008902944028897, -6.993655687767724]]), np.array([[-0.40146360986052226, -0.9003614738626398, 0.16785763712638602]])]
	}

	# FOR RCA, ASI6
	rca_asi6 = {
	'p2': [np.array([[-6.123668, 13.9552, -7.0098]]), np.array([[-.300532, -.946231, .11882]])],
	'p3': [np.array([[-6.213015, 14.166392, -7.0079]]), np.array([[-.27811, -.95650, .08807]])],
	'p4': [np.array([[-6.2683, 14.29928, -6.73526]]), np.array([[-.275515, -.961023, .0229186]])],
#	'p5': [np.array([[-5.93617, 14.29427, -6.8817]]), np.array([[-.262941, -.964801, .004606]])],
	'p5': [np.array([[-6.012922995144564, 14.32083, -6.8181226]]), np.array([[-0.268090641, -0.96323039, -0.017737533]])],
	}

	# sphere clips for rca asi4
	rca_asi4_spheres = {
	'p5': None
	}

	# FOR LAD, ASI2
	lad_asi2 = {
	'lad1': [np.array([[-12.14873847362136, 10.58989, -3.1367]]), np.array([[.63, .6467, -.4299]])], 
	'lad2': [np.array([[-12.1058, 10.605, -2.946], [-12.24216, 11.40157, -3.3217]]), np.array([[.4991, .7960, -.342366], [.63884, -.69974, .319761]])],
	'lad3': [np.array([[-12.04918, 10.51413, -2.9808], [-12.31338805, 11.38071, -3.53687]]), np.array([[.49559, .82418, -.27407], [-.337107, -.6744507, .65687]])],
	'lad4': [np.array([[-12.0736, 10.5614, -2.9405], [-12.29135, 11.42925, -3.60336]]), np.array([[.44635, .75312, -.483183], [-.2365302, -.676956, .696982]])],
	'lad5': [np.array([[-12.3028846, 10.6453, -2.9187], [-12.23315, 11.48406, -3.68404]]), np.array([[.433404, .851405, -.295415], [-.23567, -.69482, .6795]])]
	}

	# FOR LAD, ASI4
	lad_asi4 = {
	'lad1': [np.array([[-12.22223875, 10.625649302301925, -3.0766794911903426], [-12.18122189630553, 11.423373567616885, -3.575082352036909]]), np.array([[0.5219609874138667, 0.7921954572047277, -0.3162010202420827], [-0.5878579877661141, -0.40304674448105837, 0.7014102280283598]])],
	'lad2': [np.array([[-12.073312555431885, 10.62139345579502, -2.937847396851007], [-12.183116573811574, 11.554877464263935, -3.568613348863608]]), np.array([[0.5160756130120122, 0.7384903864747702, -0.43393307172721995], [-0.612329896010801, -0.42940828546696747, 0.6638227344884484]])],
	'lad3': [np.array([[-12.080526462584642, 10.586960857819236, -2.9814592093344445], [-12.201056456558787, 11.609129339800281, -3.66022598049844]]), np.array([[0.4659499051572062, 0.7897994693607993, -0.39888279491774353], [-0.5437715315845641, -0.417393972209411, 0.7280760903926368]])],
	'lad4': [np.array([[-12.017385854726623, 10.571653104861554, -2.9900184527029485], [-13.816344677148287, 12.700798817303596, -4.0591744370059235]]), np.array([[0.46793396970016604, 0.8107210244489328, -0.351808499779521], [0.9639747674317588, 0.13960718018237359, 0.22641219710168622] ])],
	'lad5': [np.array([[-12.08361770475966, 10.636157441489331, -2.9198392624436296], [-13.990413368332728, 12.553013264527257, -4.048864894266269]]), np.array([[0.45947728898589985, 0.8113527373144075, -0.3613687265918095], [0.9269190595727953, 0.2133286491762716, 0.3087263261228498]])],
	}

	# FOR LAD, ASI6
	lad_asi6 = {
	'lad1': [np.array([[-12.06955, 10.6508, -2.945514], [-13.7371712, 12.555614, -4.093051]]), np.array([[0.5519634, .775846, -.305639], [.8023378077, .400163, .442858423]])],
	'lad2': [np.array([[-11.98898, 10.57412, -3.0133259], [-13.99126, 12.571441, -3.99747] ]), np.array([[0.5644465, .65243, -.50571], [.945397, .22211, .23852]])],
	'lad3': [np.array([[-12.0323, 10.563, -3.1232], [-14.0746, 12.5727, -4.06032]]), np.array([[0.5467, .77858, -.30807], [.91761, .363155, .161603]])],
	'lad4': [np.array([[-12.0061, 10.626445, -2.958045], [-14.0991, 12.421, -4.029]]), np.array([[.563478, .756292, -.332438], [.91279, .37458, .162788]])], 
	'lad5': [np.array([[-11.9981, 10.601986, -2.9388431], [-14.1125, 12.4241255, -4.1772985]]), np.array([[.5359912, .81637837, -.2150345], [.91011496, .3682, .19006]])]
	}

	# sphere clips for lad asi4
	lad_asi4_spheres = {
	'lad1': None,
	'lad2': None,
	'lad3': None,
	'lad4': [np.array([[-11.679780025753226, 11.871200305436611, -4.34354345456884]]), np.array([0.8530250170416455])],
	'lad5': [np.array([[-11.128138767499083, 12.435758573992505, -4.456461001486442]]), np.array([1.4823622053427865])]
	}

	# sphere clips for lad asi6
	lad_asi6_spheres = {
	'lad1': [np.array([[-11.28277816, 12.2845688, -5.01866128]]), np.array([1.86445])],
	'lad2': [np.array([[-11.471234, 12.232, -4.846]]), np.array([-1.5641])],
	'lad3': [np.array([[-11.52133, 11.9257, -4.22575]]), np.array([.9483])],
	'lad4': [np.array([[-11.1808, 11.979, -4.76517]]), np.array([1.533202])],
	'lad5': [np.array([[-11.2788046, 11.9514549, -4.338]]), np.array([1.1758710])]
	}

	# encapsulate all these dicts for extensibility
	rca_total = {
	'ASI2': None,
	'ASI4': rca_asi4,
	'ASI6': rca_asi6
	}

	lad_total = {
	'ASI2': lad_asi2,
	'ASI4': lad_asi4,
	'ASI6': lad_asi6,
	}

	all_sphere_clips = None
	param_dict = None

	if 'RCA' in vessel: 
		param_dict = rca_total
		if 'ASI4' in shape: 
			all_sphere_clips = rca_asi4_spheres

	elif 'LAD' in vessel: 
		param_dict = lad_total
		if 'ASI4' in shape: 
			all_sphere_clips = lad_asi4_spheres
		elif 'ASI6' in shape: 
			all_sphere_clips = lad_asi6_spheres

	return (param_dict[shape], all_sphere_clips)


def get_reader_writer_extractor_funcs(mode): 
	if 'vtu' in mode: 
		return (return_unstructured, vtk.vtkXMLUnstructuredGridWriter(), vtk.vtkExtractGeometry())
	else: 
		return (return_polydata, vtk.vtkXMLPolyDataWriter(), vtk.vtkExtractPolyDataGeometry())


def main(): 

	# data structure is going to be dict of {suff : list of lists containing origins, normals}

	path = '/Users/alex/Documents/lab/KD-project/clipped_results_short/' 
	vessel = 'RCA/'
	shape = 'ASI4'
	mode = '.vtu'
#	mode = '.vtp'
	
	suffs = ['p1', 'p5']
#	suffs = ['lad1', 'lad2', 'lad3', 'lad4', 'lad5']

	param_dict, all_sphere_clips = get_clip_parameters(vessel, shape)

	reader, writer, extractor = get_reader_writer_extractor_funcs(mode)

	for suff, [origins, normals] in param_dict.iteritems():
		sphere_clips = None if all_sphere_clips is None else all_sphere_clips[suff]
		if suff in suffs: 
			print suff
			apply_clipping(path, vessel, shape + '/', suff,
						   normals, origins,
						   reader, writer, extractor, mode,
						   sphere_clips = sphere_clips)

		# big_operator(suff, normals, origins, sphere_clips=all_sphere_clips[suff])

if __name__ == "__main__":
	main()


