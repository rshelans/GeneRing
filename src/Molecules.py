import Chromatin
import os
import scipy
import glob


__version__="01.00.00"
__author__ ="Robert Shelansky"


class DataBase:
	def __init__(self, DATA_BASE_PATH):
		PATH_TO_DATABASE=DATA_BASE_PATH
		MOLECULE_FILES  =glob.glob(PATH_TO_DATABASE+'\\*\\*\\*\\molecule_info.txt')
		self._mols   = scipy.array([Chromatin.Molecule.info_parse(open(file)) for file in MOLECULE_FILES])
		self._name   = scipy.array([file.split(sep='\\')[-2] for file in MOLECULE_FILES]).astype(str)
		self._strain = scipy.array([file.split(sep='\\')[-4].split('_')[0] for file in MOLECULE_FILES]).astype(str) 
		self._desc   = scipy.array([set(file.split(sep='\\')[5].split('_')[1:]) for file in MOLECULE_FILES]).astype(str)


def _sample_(array, bootstraps):
	samples = (scipy.random.randint(0,len(array),size=len(array)) for i in range(bootstraps))
	for sample in samples:
		yield array[sample]

def bootstrap(array, function, bootstraps, alpha=0.95):
	strapped=scipy.array(list(map(function,_sample_(array,bootstraps))))
	sortstrapped =scipy.sort(strapped,axis=0)
	mew  = scipy.mean(sortstrapped, axis=0)
	error= scipy.array([sortstrapped[int((1-alpha)/2*bootstraps) ,:],
						sortstrapped[ bootstraps-int((1-alpha)/2*bootstraps) ,:]])
	return(mew,error,strapped)

def _rvalue_(molecules):
	scipy.sum(molecules, axis=0)/len(molecules)

def rvalue(mols):
	return(_ravlue_([mol._array_ for mol in mols]))

def _midpoint_(bubbles, length):
	return(scipy.bincount([bub.start+bub.size//2 for bub in bubbles],minlength=length))

def midpoint(bubbles, length=None, depth=1.0, gettr=None):
	assert(len(bubbles) > 0)
	if gettr == None:
		gettr=Chromatin.Molecule.getBubbles
	if isinstance(bubbles[0], Chromatin.Molecule):
		length  = len(bubbles[0])
		depth   = len(bubbles)
		bubbles = (bub for mol in bubbles for bub in gettr(mol))
	elif isinstance(bubbles[0], Chromatin.Region):
		assert(length != None)
	else:
		return(NotImplemented)	
	return(_midpoint_(bubbles,length)/depth)


def _bubble_size_(molecule):
	return( (molecule[i].size if molecule[i].isbubble else scipy.nan for i in range(len(molecule))))

def bubble_size(molecules):
	return((_bubble_size_(mol) for mol in molecules))

def test(molecules):
	return([scipy.mean([mol[i].size for mol in molecules if mol[i].isbubble]) for i  in range(len(molecules[0]))])

def smooth(array, binwidth):
	"""
	takes a list where every element is some value and preforms a sliding window average around that coordinate
	createing a new list whos len is the same as inputed.
	"""
	array         =scipy.convolve(array,scipy.ones(binwidth)/binwidth, mode='same')
	return(array)

# def rval_bubbles(bubs,factor=1):
# 	rval = scipy.zeros(len(self._mols[0]))
# 	for bub in bubs:
# 		rval[bub.start:bub.end]


def get_bubbles(molecules, start, end):
	return((bub for mol in molecules for bub in mol[start:end] if bub.isbubble))


def get_region_rval(molecules, start, end):
	return(scipy.sum(scipy.array([mol._array_[start:end] for mol in molecules]), axis=1)/(end-start))