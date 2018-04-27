'''
	parser.py

'''

import argparse
import subprocess
from pathreader import *


def gather_centerlines(model_dir):
	'''

	'''
	print '---------'
	print 'gathering centerlines from directory:'
	print model_dir
	
	# ls the target directory for *.pth 
	sp = subprocess.Popen('ls ' + model_dir + '*.pth', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	# strip the file names
	centerlineNames = [name.strip() for name in sp.stdout.readlines()]

	# additionally strip off path, .pth file ending
	names = [name[-13:-4] for name in centerlineNames]

	# user confirmation that we have thee right centerlines
	print 'gathered the following centerlines: '	
	print names

	# now that we have the centerline paths and names, read .pth files into our dictionary
	centers = {}
	for centerline, name in zip(centerlineNames, names):
		centers[name] = read_centerline(centerline)

	return (centers, names)


def parse_facenames(names, model_dir, inputFileName="SKD0050_baseline_model.vtp.facenames", ):

	file = open(model_dir+inputFileName)
	all_lines = [line.split() for line in file.readlines()]
	corresponding_faces = {}
	face_list = []
	for line in all_lines:
		if line[0] == "set":
			cur_name = line[-1][1:-1]
			if cur_name in names:
				faceID = int(line[1][19:-1])
				face_list.append(faceID) 
				corresponding_faces[cur_name] = faceID

	return (corresponding_faces, face_list)


def parse_commandline():

	# parser = argparse.ArgumentParser(description='lol')
	# parser.add_argument('--iF', action="store")
	# parser.add_argument('--kB', action="store", type=float, default=40000.0)
	# parser.add_argument('--kN', action="store", type=float, default=400.0)
	# parser.add_argument('--nbCutoff', action="store", type=float, default=0.50)
	# parser.add_argument('--m', action="store", type=float, default=12.0)
	# parser.add_argument('--dt', action="store", type=float, default=.001)
	# parser.add_argument('--n', action="store", type=int, default=1000)
	# parser.add_argument('--out', action="store", default=None ) 
	
	# args = vars(parser.parse_args())
	
	# if args['out'] == None: args['out'] = args['iF'][:-4]

	return args


if __name__ == "__main__":
	print '--------------'
	print 'testing parser.py'

	print '1) testing gathering centerlines'

	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	#print gather_centerlines(model_dir)
	centers, centerlineNames = gather_centerlines(model_dir)

	print '2) testing extraction of name-id '
	parse_facenames(centerlineNames, model_dir)
