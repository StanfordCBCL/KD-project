"""ad_utils.py


"""

# external dependencies
import numpy as np
import vtk


def produce_tagfile(total_nodes, region_ids, tagfile_name):
	"""produces a tagfile. line 0 is [total number of nodes in mesh] [1], 
	each remaining line i has 1 if node i is in the RoI, and 0 otherwise
	
	Args:
	    total_nodes (int): the total number of nodes in the complete, unclipped mesh
	    region_ids (nd_array): numpy array containing the 1-indexed node numbers contained within the RoI
	    tagfile_name (str): desired output name
	"""
	base = np.zeros(total_nodes)
	indices = region_ids - 1
	base[indices] = 1
	base = base.astype(np.int64)
	with open(tagfile_name, 'wb') as f: 
		f.write('%d\t%d\n'%(total_nodes, 1))
		base.tofile(f, sep="\n")


def main(): 
	pass

if __name__ == "__main__": 
	main()