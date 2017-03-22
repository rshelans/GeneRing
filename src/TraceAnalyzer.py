#!/usr/bin/env python
import sys
from matplotlib import pyplot as plt
import TraceView
import TraceModel
import scipy
import argparse
import glob
import os.path

__version__="01.00.00"
__author__ ="Robert Shelansky"

DEFAULT_LENGTH   =2246
DEFAULT_THRESHOLD=4
DEFAULT_SMOOTH   =10

parser = argparse.ArgumentParser(description='Analyze Chromatin Ring Trace Files.')
parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
parser.add_argument('files',
					nargs='+',
					type =str,
					help ='The list of trace files you wished to be analyzed.')

parser.add_argument('-t','--threshold',      
					default=DEFAULT_THRESHOLD,
					type   =float,
					help='The threshold distance required for two points on opposing strands to be considered as part of a linker.')

parser.add_argument('-l','--length',
					default=DEFAULT_LENGTH,
					type   =int,
					help='The expected length in basepairs for the DNA molecule.')

parser.add_argument('-s','--smooth',
					default=DEFAULT_SMOOTH,
					type   =int,
					help='The # of coordinates to include in a sliding window of the average distance between points on opposing strands.')

parser.add_argument('-u','--user',
					default=None,
					type   =str,
					help='The name of the person who completed the trace and is using this software: for record keeping.')

parser.add_argument('-o','--out_path',
					default=None,
					type   =str,
					help='The name of the folder you wish to save the output to.')

parser.add_argument('-p', '--plotless',
					action='store_true')

parser.add_argument('-i','--imgres',
					default = scipy.NAN,
					type    = scipy.float16,
					help    ='The image resolution of the raw image used to make the trace.')
parser.add_argument('-d','--directory',
					action='store_false',
					default='true')


args      = parser.parse_args()
args.files=[f for g in args.files for f in glob.glob(g)]
if len(args.files) == 0:
	sys.stderr.write('No Trace Files Found. Check the file path.')
	sys.exit()


params    = vars(args)
params['version']=__version__
params['title']=os.path.basename(sys.argv[0]).split('.')[0]


model = TraceModel.TraceModel(**params)
if args.plotless:
	for i in range(len(args.files)):
		model.seek(i)
		model.find_midpoint()
		model.analyze()
		if model.molecule['molecule'] is None:
			print(i, model.context['path'])
		else:
			model.save()
else:
	view  = TraceView.TraceView (plt, model)
	view.show()