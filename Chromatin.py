import scipy
import collections
import bisect
import sys

__version__="01.00.00"
__author__ ="Robert Shelansky"

BUBBLE = int(True)
LINKER = int(False)

ZRS =(0,78)
UAS1=(281,303)
UAS2=(389,411)
TATA=(551,557)
TSS = 607
ORF =(652,2055)
NUCLEOSOME_SIZE=147
NUCLEOSOME_CUTOFF=90

#Robert's Regions
N1REGION=[UAS2[1]+24,TSS+NUCLEOSOME_SIZE//2+5]
N3REGION=[ZRS[0],UAS2[0]-35]
N2REGION=[UAS2[0]-83,UAS2[1]+NUCLEOSOME_CUTOFF+10]
#Saketh's Regions
sN1REGION=[490,580]
sN2REGION=[UAS2[0]-83,UAS2[1]+NUCLEOSOME_CUTOFF+10]
sN3REGION=[ZRS[0],UAS2[0]-35]


class Region:
		"""
		Class Molecule represents individual bubbles or linkers
		"""

		def __init__(self, bubble, start, end, size):
			self.isbubble= int(bubble)
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
		self._array_  = scipy.array(list("".join([str(reg.isbubble)*reg.size for reg in self]))).astype(int)

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




def get_configuration(molecule):
	"""
	get_configuration takes a molecule object and returns its estimation of the 
	configuration of the molecuels promoter as defined by 
	
	'Nucleosomal promoter variation generates gene expression noise'
	Christopher R. Brown and Hinrich Boeger1

	The return type is the configuration Number.
	"""

	#Definitions for molecule Regions
	n1   =    [bub for bub in molecule.getOverlap(N1REGION[0],N1REGION[1],NUCLEOSOME_CUTOFF) if bub.isbubble]
	n2   =    [bub for bub in molecule.getExc(UAS2[0],UAS2[1])  if (bub.isbubble and bub.size >NUCLEOSOME_CUTOFF) ]
	n3   =    [bub for bub in molecule.getOverlap(N3REGION[0],N3REGION[1],NUCLEOSOME_CUTOFF) if bub.isbubble]
	return(Configuration(n1,n2,n3))
