'''
	parser.py

'''

import argparse
import subprocess


def gather_centerlines(path_direct):
	sp = subprocess.Popen('ls ' + path_direct + '*.pth', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	paths = [name.strip() for name in sp.stdout.readlines()]
	return (len(paths), paths)



if __name__ == "__main__":
	path_direct = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	print gather_centerlines(path_direct)
