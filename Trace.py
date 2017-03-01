import scipy
import scipy.spatial
import itertools
import Chromatin

__version__="01.00.00"
__author__ ="Robert Shelansky"

class Trace:
	"""
	Trace represents an individual trace of a molecule from a image file. It can be instantiated by giving it
	a reference to a list of tupples of coordinates from starting at the end of the molecule (the fork) and leading
	to the beginning of the molecule. It contains all of the data required to do analysis of an individual trace.
	as well as helper functions which both do the analysis and manipulate an individual "realization of the fit" -- 
	this simply meaning for each coordinate whether its part of a linker or not.
	"""
	def __init__(self, trace):
		"""
		Takes a list of coordinate tupples and computes metrics required for realizing a specific bubble linker path.
		usable metrics are as follows.
		_trace:
			#array of x,y coordinates of on single _trace
		_ld:
			#distance between succesive points linked diff (ld) and the distance be all points as a matrix (d)
			#index 0 refers to distance between 0,1 in _trace, distance
			#index -1 refers to distance between -2,-1 in _trace
		_cld:
			#cumulative distance between coordinates starting at 0,1
			#there is no index 0 
			#index i refers to the distance traveled to get to index i+1 in _trace	
		_ll:
			#length of the whole molecule in the coordinate system
		_d:
			#distance between every point and every other
		"""
		self._trace      =scipy.array(trace)
		self._ld         =scipy.array([scipy.spatial.distance.euclidean(i,j) for i,j in zip(self._trace[:-1],self._trace[1:])])
		self._cld        =scipy.concatenate(([0],scipy.cumsum(self._ld)))
		self._ll         =scipy.sum(self._ld)
		self._d          =scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(self._trace,'euclidean'))
		self._d          =scipy.ma.masked_where(self._d==0,self._d) ##mask self distances

	def __getitem__(self, arg):
		"""
		returns the coordinate at index arg. (   [x,y]  )
		"""
		return (self._trace[arg])

	def __len__(self):
		"""
		returns the length in coordinates of the trace. The number of coordinate entries.
		"""
		return(len(self._trace))

	
	def partition(self, midpoint):
		"""
		returns two lists forward and reverse who have length len(self)=len(forward) +len(reverse)
		each index contains the distance of the closest coordinate of the other strand of the molecule as defined by midpoint
		"""
		forward=scipy.amin(self._d[0:midpoint+1,midpoint+1:len(self)],axis=1)
		reverse=scipy.amin(self._d[0:midpoint+1,midpoint+1:len(self)],axis=0)
		return forward, reverse

	def smooth(self, mol, smooth):
		"""
		takes a list where every element is some value and preforms a sliding window average around that coordinate
		createing a new list whos len is the same as inputed.
		"""
		mol         =scipy.convolve(mol,scipy.ones(smooth)/smooth,mode='same')
		return(mol)

	def segment(self, mol, threshold, end):
		"""
		takes a list whos values are some score of whether that specific coordinate is a linker.
		this is typically used on the lists returned from partition but does not have to be.
		All streches of elements whos threshold is above or below some value are then determined.
		and a list of regions and there corresponding label are returned. It also ignores the last "end"
		coordinates of mol.

		regeions = [(0,end),(start,end),...,(start,end),(start,len(mol))]
		labels   = [1,0,...1,0]  1 for greater then threshold 0 for below.

		"""
		mask        =mol > threshold
		breaks      =scipy.concatenate(( [0],scipy.argwhere(scipy.ediff1d(mask[:-end]))[:,0]+1, [len(mol)] ))
		regions     =scipy.array([(s,e) for s,e in zip(breaks[:-1],breaks[1:])])
		labels      =scipy.array([int(mask[s]) for s,e in regions])
		return(regions,labels)

	def zip(self, fregions, flabels, rregions, rlabels):
		"""
		takes the definitions of forward regions and reverse regions and returns the concatenate version.
		This is simply the region definitions for the molecule in its entiery
		"""
		regions = scipy.concatenate((fregions ,rregions))
		labels  = scipy.concatenate((flabels,rlabels))
		return(regions,labels)

	def msd_of_linker_end_points(self,fregions,flabels,rregions,rlabels):
		"""
		takes definitions of regions and labels for forward and reverse
		strands assuming each region in the forward direction coorosponds to a region
		in the reverse direction and returns the mean squaared distance
		of the start and end points of each region. 
		"""
		if len(fregions) != len(rregions):
			return(float('Inf'))
		flinks=fregions[scipy.logical_not(flabels)]
		rlinks=rregions[scipy.logical_not(rlabels)]
		s=sum([ self._d[t1,t2] for t1,t2 in zip(flinks[:,0]  ,rlinks[:,1] -1)  if self._d[t1,t2]])   /len(flinks)
		f=sum([ self._d[t1,t2] for t1,t2 in zip(flinks[:,1]-1,rlinks[:,0])     if self._d[t1,t2]])   /len(flinks)
		return((s+f)/2)

	def msd_of_region_sizes(self,fregions,flabels,rregions,rlabels):
		"""
		takes definitions of regions and labels for forward and reverse
		strands assuming each region in the forward direction coorosponds to a region
		in the reverse direction and returns the mean squared distance between the sizes of 
		each region.
		"""
		if len(fregions) != len(rregions):
			return(float('Inf'))
		flen = self._cld[fregions[-1,-1]-1] - self._cld[fregions[0,0]]
		rlen = self._cld[rregions[0,-1] -1] - self._cld[rregions[-1,0]]
		dif  = sum([((self._cld[ff-1]-self._cld[fs])/flen-(self._cld[rf-1]-self._cld[rs])/rlen) **2            for (fs,ff),(rs,rf) in  zip(fregions,rregions)])
		dif  = dif*(flen+rlen)/2
		return(dif)

	def sd_segment_size(self,fregions,flabels,rregions,rlabels):
		"""
		takes definitions of regions and labels for forward and reverse
		strands assuming each region in the forward direction coorosponds to a region
		in the reverse direction and returns how similar in size each fragment is.
		"""
		if len(fregions) != len(rregions):
			return(float('Inf'))
		flen = self._cld[fregions[-1,-1]-1] - self._cld[fregions[0,0]]
		rlen = self._cld[rregions[0,-1] -1] - self._cld[rregions[-1,0]]
		return((flen - rlen)**2)

	def regionify(self, midpoint, threshold=4, smooth=10, end=5):
		"""
		returns a pair of tuples that is a best guess at regions given a specific
		midpoint threshold end and smooth factor. A convinience function which calls other
		trace methods.
		"""
		a,b  = self.partition( midpoint)
		f,fl = self.segment  (self.smooth(a      , smooth), threshold=threshold, end= end)
		r,rl = self.segment  (self.smooth(b[::-1], smooth), threshold=threshold, end= end)
		r    = len(self)-r[:,[1,0]]
		rl   = rl
		return((f,fl),(r,rl))


	def _midpoint(self, left, right):
		"""
		given two coordinates returns the middle coordinate in terms of actual distance.
		"""
		return(scipy.searchsorted(self._cld, self._cld[left]+(self._cld[right]-self._cld[left]) / 2 ))

	# def find_midpoint(self, guess=None ,threshold=4, sensitivity=10, end=5):
	# 	guess           = guess or scipy.searchsorted(self._cld, (self._cld[-1]/2))
	# 	(ff,fl),(rf,rl) = regionify(self, guess, threshold= threshold, sensitivity=sensitivity, end=end)
	# 	i=min(len(ff),len(rf))-2
	# 	guess = _midpoint(self, ff[i][1],rf[i][0])
	# 	return(guess)
	def midpoint(self,guess=None,threshold=4, smooth=10, end=5):
		"""
		Takes some perameters and returns the best guess of the midpoint.
		It does this iteritivly. It attemps to define regions based on midpoint. Then looks at the
		second to last set of regions in the shorter strand and attempts to recalculate the midpoint 
		given that region and its coorrosponding region in the longer strand. does this until convergence.
		or (until it attempts number of regions found in one long molecule)--uppper bound not super important.
		helps in absence of convergence.
		"""
		guess           = guess or scipy.searchsorted(self._cld, (self._cld[-1]/2))
		(ff,fl),(rf,rl) = self.regionify(guess, threshold= threshold, smooth=smooth, end=end)
		for i in range(min(len(ff),len(rf))):
			i         = min(len(ff),len(rf))-2
			new_guess = self._midpoint(ff[i][1],rf[i][0])
			if guess==new_guess:
				return(guess)
			guess = new_guess
			(ff,fl),(rf,rl) = self.regionify(guess, threshold= threshold, smooth=smooth, end=end)
		return(guess)

	def edgebuffer(self, threshold, smooth):
		"""
		Calculates how many coordinates to ignore on the end by determining
		the ceiling of the minimum number of coordinates to meet threshold
		"""
		return(int(scipy.ceil(threshold/min(self._ld))))

	def solve_molecule(self, midpoint, threshold=4, smooth=10, end=5):
		"""
		given specific settings produces a list of objects which represent a realization of a trace.
		It is again a convinecnce funciton like regionify 
		"""
		molecule       =self.smooth    (scipy.concatenate(self.partition(midpoint)),smooth=smooth)
		(fr,fl),(rr,rl)=self.regionify (midpoint, threshold=threshold , smooth=smooth,end=end)
		regions, labels=self.zip       (fr,fl,rr[::-1],rl[::-1])
		return (midpoint,molecule,(fr,fl),(rr,rl),(regions,labels))

	def label(self,fl,rl):
		"""
		given two lists of regions or labels: it returns a list of of length len(trace) whos
		values are which region it belongs to if you attempt to zip of the molecule from the end.
		"""
		return(scipy.concatenate((list(range(len(fl))),list(reversed(range(len(rl)))))))

	def scale(self,length):
		"""
		calculates the number of basepairs per coordinate distance.
		"""
		return(length/self._cld[-1])

	def moleculify(self,fr,fl,rr,rl,length):
		"""
		takes a representation of a trace fir region definitions and labels and a length in basepairs
		of a molecule and returns a Chromatin.molecule version.
		"""
		#mean length
		if len(fl)!=len(rl) or (sum((fl + rl)==1) > 0):
			return(None)
		
		region_lengths       = scipy.array([sum((self._cld[r1[1]-2] - self._cld[r1[0]], self._cld[r2[1]-2] - self._cld[r2[0]]))/2  for r1,r2 in zip(fr,rr)])
		exclusive_end_pts    = scipy.ceil(length * scipy.ndarray.round(scipy.cumsum(region_lengths)/sum(region_lengths),decimals=3))

		inclusive_start_pts  = scipy.concatenate(([0],exclusive_end_pts[:-1]))
		
		regions  = scipy.array([(s,e) for s,e in zip(inclusive_start_pts,exclusive_end_pts)])
		molecule=Chromatin.Molecule([ Chromatin.Region(l,length-e,length-s,e-s) for (s,e),l in reversed(list(zip(regions,fl)))])
		return(molecule)
