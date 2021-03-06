import scipy
import collections
import bisect
import sys


__version__="01.00.00"
__author__ ="Robert Shelansky"

BUBBLE = int(True)
LINKER = int(False)

class Region:
		"""
		Class Molecule represents individual bubbles or linkers
		"""

		def __init__(self, bubble, start, end, size):
			self.isbubble= bool(bubble)
			self.start = int(start)
			self.end   = int(end)
			self.size  = int(size)
			assert(self.end-self.start == self.size)

		def __str__(self):
			return(",".join([str(self.isbubble),str(self.start),str(self.end),str(self.size)]))

		def __len__(self):
			return(self.size)

class Molecule:
	"""
	Class Molecule is a representation of Individual Traces of chromatin Molcules off the EM.
	One molecule can be indexed for region objects these are objects which represent bubbled.
	And not bubbled DNA.
	i.e.

	mol[10]    # one linker or bubble which ocupies that space
	mol[10:20] #returns all linkers and bubbles which occupy that space

	-Assumes the whole molecule is made of regions.
	-includes both bubbles and linkers for now.
	-stops are inclusive
	"""
	def __init__(self,region_list):
		self.regions  = sorted(region_list,key=lambda x:x.start)
		self._ends_   = [region.end for region in self.regions]
		self._starts_ = [region.start for region in self.regions]
		self._array_  = Molecule.array(self.regions,len(self)).astype(int)

	def getBubbles(self):
		return([region for region in self.regions if region.isbubble])

	def getLinkers(self):
		return([region for region in self.regions if not region.isbubble])
		
	def __getitem__(self, arg):
		"""
		option (1) Returns All regions within the slice if step # not present
		           (exception: has to be a slice obj with no third arg)

		option (2) Regions that contain all the bases called on will be returned 
		           (exception: a third integer arg (1) must be passed in)

		option (3) Region with the given base will be returned
		           (exception: only one integer argument should be passed in)
		"""	
		if isinstance(arg, slice):	
			return (self.getExc(arg.start,arg.stop))
		else:
			return (self.getReg(arg))

	def getReg(self, pos):
		"""getReg(pos) returns the region that occpies position"""
		return(self.regions[bisect.bisect_left(self._ends_,pos)])


	def getExc(self, beg=None,fin=None):
		"""
		 getExc(beg,fin) returns all regions that has a part of itself (or entirely) in the slice

		"""
		start = bisect.bisect_left(self._ends_,beg) if beg else 0
		stop = bisect.bisect_left(self._ends_,fin) if fin else self._ends_[-1]
		return(self.regions[start:stop+1])

	def getInc(self, beg=None,fin=None):
		"""
		getInc(beg, fin) returns only the region that has all it's bases included in the slice
		
		"""
		
		start = bisect.bisect_right(self._starts_,beg) if beg else 0
		stop = bisect.bisect_left(self._ends_, fin) if fin else self._ends_[-1]
		return (self.regions[start:stop])

	def getOverlap(self,beg=None,fin=None,overlap=90):
		"""getOverlap returns only regions which contain x number of base overlap with the slice"""
		return ([region for region in self.getExc(beg,fin) if (min(fin,region.end)-max(beg,region.start)) > overlap]) 

	def __len__(self):
		#end - Beginning (0)
		return(self._ends_[-1])

	def __iter__(self):
		for region in self.regions: 
			yield region

	def write(self, file=sys.stdout):
		"""
		Writes a Molecule to A file in the following format:
		## - comments
		isbubble\tstart\tend
		0\t10\t100
		...
		...
		...
		"""
		for region in self:
			print("{}\t{}\t{}".format(region.isbubble,region.start,region.end),file=file) 


	@staticmethod
	def array(regions, length):
		"""
		"""
		array   = scipy.empty(length)
		array[:]= scipy.NAN
		for reg in regions:
			array[reg.start:reg.end] = reg.isbubble
		return (array)

	@staticmethod
	def filter(molecule, filter_func):
		"""
		Beta version:
		Uses the filter_func which must return a boolean when being passed an individual region to filter
		out certain regions. !!!DO NOT FILTER LINKERS AND BUBBLES AT THE SAME TIME.!!! All positions that were removed
		are then interpolated by the surrounding values.
		"""
		filtered       = [Region(r.isbubble,r.start,r.end,len(r)) for r in molecule if not filter_func(r)]
		filtered       = Molecule.array(filtered,len(molecule))
		#Selects all removed bases and attempts to interpolate their value from surrounding bases.
		#http://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
		nans           = scipy.isnan(filtered)
		f              = lambda x: x.nonzero()[0]
		filtered[nans] = scipy.interp(f(nans),f(~nans),filtered[~nans])
		#After interpolation the regions are reconstructed an a new molecule is created with specific regions removed.
		pivots         = scipy.argwhere(scipy.diff(filtered))[:,0]+1
		pivots         = scipy.concatenate(([0],pivots,[len(molecule)]))
		new_regions    = [Region(bool(filtered[start]),start,stop,stop-start) for start,stop in zip(pivots[:-1],pivots[1:])]
		return(Molecule(new_regions))

	@classmethod
	def read(cls,file):
		"""
		takes a .mol file and buils a molecule object
		"""
		data   =scipy.genfromtxt(file,comments='#').astype(int)
		regions=[Region(i,s,e,e-s) for i,s,e in data]
		return(Molecule (regions) )
	
	@classmethod
	def info_parse(cls,  file):
		"""
		Depricated: This is an old file format Prefered is to use the new file
		Format.
		input is each molecule_info.txt file: 
			#base_in_one_strand,base_in_other,bubble_pos,linker_pos
		output is a list 
			#bubble/linker,pos_in_molecule,start_base,end_base,size
		Could be written infinitly btter why worry about it now
		"""
		features=collections.defaultdict(lambda:[])
		bubble_list=[]
		for i,line in enumerate(file):
			line=line.strip().split()
			if int(line[2]) != 0:		
				features["BUBBLE"+line[2]].append(line + [i])
			else:
				features["NOTBUBBLE"+line[3]].append(line +[i])
		
		for feat in features.keys():
			is_bubble=int(feat.startswith("BUBBLE"))
			feat = [x[4] for x in features[feat]]
			start= feat[0]
			end  = feat[len(feat)-1]
			size = end-start+1	
			bubble_list.append(Region(is_bubble,start,end+1,size))
		return Molecule(bubble_list)




