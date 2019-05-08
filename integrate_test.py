"""integrate_test.py



"""
import numpy as np
import vtk
from vtk.util import numpy_support as nps	

from AneurysmGeneration.utils.batch import return_polydata
from AneurysmGeneration.utils.slice import extract_points


def batched_triangle_areas(points1, points2, points3): 
	"""compute the surface area for a bunch of triangles
	
	Args:
	    points1 (ndarray): shape (NoC, 3), represents 1st point of each triangle
	    points2 (ndarray): shape (NoC, 3), represents 2nd point of each triangle
	    points3 (ndarray): shape (NoC, 3), represents 3rd point of each triangle
	
	Returns:
	    S: ndarray of shape (NoC, 1), representing area of each triangle
	"""
	u = points2 - points1
	v = points3 - points1

	S = .5*np.linalg.norm(np.cross(u, v), axis=-1)

	return S


def main(): 
	file_path = "clipped_results_short/RCA/ASI2/p1.vtp"

	polydata = return_polydata(file_path)

	NoP = polydata.GetNumberOfPoints()
	NoC = polydata.GetNumberOfCells()

	print(NoC)
	print('----')
	raw_vtawss = nps.vtk_to_numpy(polydata.GetPointData().GetArray('vTAWSS'))
	_, raw_points = extract_points(polydata)

	local_averages = np.zeros(NoC)
	points1, points2, points3 = np.zeros((NoC, 3)), np.zeros((NoC, 3)), np.zeros((NoC, 3))

	for i in range(NoC):
		
		cell_pt_ids = [int(polydata.GetCell(i).GetPointIds().GetId(j)) for j in range(3)]
		local_averages[i] = np.mean([raw_vtawss[pt_id] for pt_id in cell_pt_ids])
		for arr, pt_id in zip([points1, points2, points3], cell_pt_ids): 
# 			print(pt_id)
			arr[i] = raw_points[pt_id]

	cell_areas = batched_triangle_areas(points1, points2, points3)
	weighted_average = np.average(local_averages, weights=cell_areas)

	print(np.mean(raw_vtawss))
	print(np.mean(local_averages))
	print(weighted_average)
	print(np.sum(cell_areas))


if __name__ == "__main__": 
	main()