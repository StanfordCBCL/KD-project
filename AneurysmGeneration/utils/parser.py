'''
	parser.py

'''

import argparse
import subprocess


def gather_centerline_names(path_direct):
	sp = subprocess.Popen('ls ' + path_direct + '*.pth', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	paths = [name.strip() for name in sp.stdout.readlines()]
	return paths


def gather_centerlines(model_dir):
	centerlineNames = gather_centerline_names(model_dir)
	centers = []
	print model_dir
	print centerlineNames
	for centerline in centerlineNames:
		centers.append(read_centerline(centerline))

def parse_commandline():

	parser = argparse.ArgumentParser(description='lol')
	parser.add_argument('--iF', action="store")
	parser.add_argument('--kB', action="store", type=float, default=40000.0)
	parser.add_argument('--kN', action="store", type=float, default=400.0)
	parser.add_argument('--nbCutoff', action="store", type=float, default=0.50)
	parser.add_argument('--m', action="store", type=float, default=12.0)
	parser.add_argument('--dt', action="store", type=float, default=.001)
	parser.add_argument('--n', action="store", type=int, default=1000)
	parser.add_argument('--out', action="store", default=None ) 
	
	args = vars(parser.parse_args())
	
	if args['out'] == None: args['out'] = args['iF'][:-4]

	return args



	return np.array(centers)
if __name__ == "__main__":
	path_direct = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	print gather_centerline_names(path_direct)
