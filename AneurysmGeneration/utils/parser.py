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

	# user confirmation that we have thee right centerlines
	print 'gathered the following centerlines: '
	names = [name[-13:-4] for name in centerlineNames]
	print names

	centers = []
	for centerline in centerlineNames:
		centers.append(read_centerline(centerline))

	return (centers, names)


def parse_facenames(names, model_dir, inputFileName="SKD0050_baseline_model.vtp.facenames", ):

	file = open(model_dir+inputFileName)
	all_lines = [line.split() for line in file.readlines()]
	corresponding_faces = []
	for line in all_lines:
		if line[0] == "set" and (line[-1][1:-1] in names):
			corresponding_faces.append(int(line[1][19:-1])) 

	return corresponding_faces


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


if __name__ == "__main__":
	print '--------------'
	print 'testing parser.py'

	print '1) testing gathering centerlines'

	model_dir = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	#print gather_centerlines(model_dir)
	centers, centerlineNames = gather_centerlines(model_dir)

	print '2) testing extraction of name-id '
	parse_facenames(centerlineNames, model_dir)
