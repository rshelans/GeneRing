import Chromatin
import os
import scipy
import glob
import matplotlib
import seaborn
__version__="01.00.00"
__author__ ="Robert Shelansky"
class Analysis:
	def __init__(self, files):
		self.__add_mols__(files)

	def __add_mols__(self, files):
		self._mols  = scipy.array([Chromatin.Molecule.info_parse(open(file)) for file in files ])
		self._name  = scipy.array([file.split(sep='\\')[-2] for file in files]).astype(str)
		self._strain= scipy.array([file.split(sep='\\')[-4].split('_')[0] for file in files]).astype(str) 
		self._desc  = scipy.array([set(file.split(sep='\\')[5].split('_')[1:]) for file in files]).astype(str)

class Subset:
	def __init__(self,mols,color='red',label='label'):
		self._label  = label
		self._color  = color
		self._mols   = mols
		self._bubs   = scipy.array([bub for mol in self._mols for bub in mol if bub.isbubble])
		self._asarray= scipy.array([mol._array_ for mol in mols]) 
		self._rval   = scipy.sum(self._asarray,0) / len(self._asarray)
		self._config       = scipy.array([Chromatin.get_configuration(mol) for mol in self._mols])
		self._config_ints  = scipy.array([config.configuration for config in self._config]).astype(int)
		self._config_hist   = scipy.bincount(self._config_ints)/len(self._config_ints)
		self._config_order  = scipy.array([i for  i in sorted(range(len(self._mols)),   key=lambda x:self._config_ints[x])])
		
		#n-3,n-2,n-1, Occupancy
		self._nuc_hist      = scipy.array((sum(self._config_hist[[0,2,6,1]]),sum(self._config_hist[[0,1,3,5]]),sum(self._config_hist[[0,2,3,4]]) )  )

		#Create a array of the moleucle rvalue only with nuc size things		
		self._nucsarray      = scipy.array(list("".join([str(int(reg.size>=90)*reg.isbubble)*reg.size  for mol in self._mols for reg in mol]))).astype(int)
		self._nucsarray.shape= [len(self._mols),len(self._mols[0])]
		self._nucs_rval      =scipy.sum(self._nucsarray,0)/len(self._nucsarray)
		
		#Create a array of the moleucle rvalue only with prenuc size things		
		self._prenucsarray      = scipy.array(list("".join([str(int(reg.size<90)*reg.isbubble)*reg.size  for mol in self._mols for reg in mol]))).astype(int)
		self._prenucsarray.shape= [len(self._mols),len(self._mols[0])]
		self._prenucs_rval      = scipy.sum(self._prenucsarray,0)/len(self._prenucsarray)

		#a boolean array to select mols that have specific configs
		self._s    ={i:scipy.array([bool(Chromatin.get_configuration(mol).configuration == i) for mol in self._mols]) for i in range(8) }
		self._n1   =scipy.amax([self._s[i] for i in [0,2,3,4]],0)
		self._n2   =scipy.amax([self._s[i] for i in [0,1,3,5]],0)
		self._n3   =scipy.amax([self._s[i] for i in [0,2,6,1]],0)
		self._n1_n2=scipy.amax([self._s[i] for i in [0,3]],0)
		self._n1_n3=scipy.amax([self._s[i] for i in [0,2]],0)
		self._n2_n3=scipy.amax([self._s[i] for i in [0,1]],0)

		self._n1_not_n2=scipy.amax([self._s[i] for i in [2,4]],0)
		self._n1_not_n3=scipy.amax([self._s[i] for i in [3,4]],0)
		self._n2_not_n3=scipy.amax([self._s[i] for i in [3,5]],0)
		self._n2_not_n1=scipy.amax([self._s[i] for i in [1,5]],0)
		self._n3_not_n1=scipy.amax([self._s[i] for i in [1,6]],0)
		self._n3_not_n2=scipy.amax([self._s[i] for i in [2,6]],0)

		self._n1_n2_n3        =scipy.amax([self._s[i] for i in [0]],0)
		self._n3_not_n1_not_n2=scipy.amax([self._s[i] for i in [6]],0)
		self._n2_not_n1_not_n3=scipy.amax([self._s[i] for i in [5]],0)
		self._n1_not_n2_not_n3=scipy.amax([self._s[i] for i in [4]],0)

		self._condusive   = scipy.amax([self._s[i] for i in [2,4,6,7]],0)
		self._incondusive = scipy.amax([self._s[i] for i in [0,1,3,5]],0)

	def rval(self,mols):
		return(scipy.sum([mol._array_ for mol in mols],0)/len(mols))

	def rval_smooth(self,mols,binwidth):
		rval = scipy.sum([mol._array_ for mol in mols],0)/len(mols)
		return(scipy.array([sum(rval[i-binwidth//2:i+binwidth//2])/binwidth for i in range(0,len(self._mols[0])-binwidth)]))

	def midpoint(self,bubbles,factor=1):
		return(scipy.bincount([bub.start+bub.size//2 for bub in bubbles])/factor)

	def midpoint_smooth(self,bubbles,binwidth,factor=1):
		midpoints = scipy.bincount([bub.start+bub.size//2 for bub in bubbles])/factor
		return(scipy.array([sum(midpoints[i-binwidth//2:i+binwidth//2])/binwidth for i in range(0,len(self._mols[0])-binwidth)]))


	def smooth(self, array, binwidth):
		return(scipy.array([sum(array[i-binwidth//2:i+binwidth//2]) for i in range(0,len(array)-binwidth)]))

	def rval_bubbles(self,bubs,factor=1):
		rval = scipy.zeros(len(self._mols[0]))
		for bub in bubs:
			rval[bub.start:bub.end]
		"""
		#dependensies of Nucleosome position on configuration
		self._n2_on_n1_pos = scipy.zeros(len(self._asarray[0]))
		self._n2_on_n1_bubble =  [reg for config in self._n2 for reg in config.n1]
		for reg in self._n2_on_n1_bubble:
				self._n2_on_n1_pos[reg.start:reg.end] +=1
		
		self._n2_on_n1_pos = self._n2_on_n1_pos/len(self._n2_on_n1_bubble)

		self._n2_on_n3_pos = scipy.zeros(len(self._asarray[0]))
		self._n2_on_n3_bubble =  [reg for config in self._n2 for reg in config.n3]
		for reg in self._n2_on_n3_bubble:
				self._n2_on_n3_pos[reg.start:reg.end] +=1
		
		self._n2_on_n3_pos = self._n2_on_n3_pos/len(self._n2_on_n3_bubble)

		self._n1_on_n3_pos = scipy.zeros(len(self._asarray[0]))
		self._n1_on_n3_bubble =  [reg for config in self._n1 for reg in config.n3]
		for reg in self._n1_on_n3_bubble:
				self._n1_on_n3_pos[reg.start:reg.end] +=1
		
		self._n1_on_n3_pos = self._n1_on_n3_pos/len(self._n1_on_n3_bubble)

		self._n3_on_n1_pos = scipy.zeros(len(self._asarray[0]))
		self._n3_on_n1_bubble =  [reg for config in self._n3 for reg in config.n1]
		for reg in self._n3_on_n1_bubble:
				self._n3_on_n1_pos[reg.start:reg.end] +=1
		
		self._n3_on_n1_pos = self._n3_on_n1_pos/len(self._n3_on_n1_bubble)

		self._n3_on_n2_pos = scipy.zeros(len(self._asarray[0]))
		self._n3_on_n2_bubble =  [reg for config in self._n3 for reg in config.n2]
		for reg in self._n3_on_n2_bubble:
				self._n3_on_n2_pos[reg.start:reg.end] +=1
		
		self._n3_on_n2_pos = self._n3_on_n2_pos/len(self._n3_on_n2_bubble)

		#every config that has a nuc at a specific position
		self._not_n1 =  [config  for c,config in zip(self._config_ints,self._config) if int(c) not in [0,2,3,4]]
		self._not_n2 =  [config  for c,config in zip(self._config_ints,self._config) if int(c) not in [0,1,3,5]]
		self._not_n3 =  [config  for c,config in zip(self._config_ints,self._config) if int(c) not in [0,2,6,1]]


		#dependensies of Nucleosome position on configuration
		self._not_n2_on_n1_pos = scipy.zeros(len(self._asarray[0]))
		self._not_n2_on_n1_bubble =  [reg for config in self._not_n2 for reg in config.n1]
		for reg in self._not_n2_on_n1_bubble:
				self._not_n2_on_n1_pos[reg.start:reg.end] +=1
		
		self._not_n2_on_n1_pos = self._not_n2_on_n1_pos/len(self._not_n2_on_n1_bubble)

		self._not_n2_on_n3_pos = scipy.zeros(len(self._asarray[0]))
		self._not_n2_on_n3_bubble =  [reg for config in self._not_n2 for reg in config.n3]
		for reg in self._not_n2_on_n3_bubble:
				self._not_n2_on_n3_pos[reg.start:reg.end] +=1
		
		self._not_n2_on_n3_pos = self._not_n2_on_n3_pos/len(self._not_n2_on_n3_bubble)

		self._not_n1_on_n3_pos = scipy.zeros(len(self._asarray[0]))
		self._not_n1_on_n3_bubble =  [reg for config in self._not_n1 for reg in config.n3]
		for reg in self._not_n1_on_n3_bubble:
				self._not_n1_on_n3_pos[reg.start:reg.end] +=1
		
		self._not_n1_on_n3_pos = self._not_n1_on_n3_pos/len(self._not_n1_on_n3_bubble)

		self._not_n3_on_n1_pos = scipy.zeros(len(self._asarray[0]))
		self._not_n3_on_n1_bubble =  [reg for config in self._not_n3 for reg in config.n1]
		for reg in self._not_n3_on_n1_bubble:
				self._not_n3_on_n1_pos[reg.start:reg.end] +=1
		
		self._not_n3_on_n1_pos = self._not_n3_on_n1_pos/len(self._not_n3_on_n1_bubble)

		self._not_n3_on_n2_pos = scipy.zeros(len(self._asarray[0]))
		self._not_n3_on_n2_bubble =  [reg for config in self._not_n3 for reg in config.n2]
		for reg in self._not_n3_on_n2_bubble:
				self._not_n3_on_n2_pos[reg.start:reg.end] +=1
		
		self._not_n3_on_n2_pos = self._not_n3_on_n2_pos/len(self._not_n3_on_n2_bubble)

		#every region that is considered a specific bubble
		self._n1_bubble = [reg for config in self._config for reg in config.n1]
		self._n2_bubble = [reg for config in self._config for reg in config.n2]
		self._n3_bubble = [reg for config in self._config for reg in config.n3]

		#an array indicating the probability of a nuc lieing in a specific region
		self._n1_pos = scipy.zeros(len(self._asarray[0]))
		for reg in self._n1_bubble:
			self._n1_pos[reg.start:reg.end] +=1
		
		self._n2_pos = scipy.zeros(len(self._asarray[0]))
		for reg in self._n2_bubble:
			self._n2_pos[reg.start:reg.end] +=1
		
		self._n3_pos = scipy.zeros(len(self._asarray[0]))
		for reg in self._n3_bubble:
			self._n3_pos[reg.start:reg.end] +=1

		self._n1_pos = self._n1_pos/len(self._n1_bubble)
		self._n2_pos = self._n2_pos/len(self._n2_bubble)
		self._n3_pos = self._n3_pos/len(self._n3_bubble)

		#Find a better way of looking at positioning maybe midpoints?
		self._mids_all = scipy.array([bub.start+bub.size//2 for bub in self._bubs])
		self._mids_nucs= scipy.array([bub.start+bub.size//2 for bub in self._bubs if bub.size > 90 and bub.size < 204])

"""




def build_config_plot(plt):
	fig = plt.figure()
	ax1 = plt.subplot2grid((10,4), (0,0), colspan=4,rowspan=8)
	ax2=  plt.subplot2grid((10,4), (8 ,0),colspan=4,rowspan=2,sharex=ax1)
	ax1.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
	ax2.tick_params(axis='y', which='both',left='off',right='on',labelright='on', labelleft='off')
	circles=[[matplotlib.patches.Circle((i,3-n),0.20,color='black') for n,c in enumerate(config) if c] for i,config in Chromatin.Configuration.MAP.items()]
	circles=[c for l in circles for c in l]
	boxes  =[matplotlib.patches.Rectangle((i -0.75/2 ,0), 0.75, 4, lw=4,edgecolor='black',fill=False) for i in range(8)]
	[ax2.add_artist(artist) for artist in boxes]
	[ax2.add_artist(artist) for artist in circles]
	ax2.set_yticklabels(["5'",'','N -3','','N -2','','N -1','',"3'"])

	ax2.set_ylim([0,4])
	ax2.set_xlim([-0.5,7.5])
	return(fig,(ax1,ax2))

	
def build_nuc_config_plot(plt):	
	fig = plt.figure()
	ax1 = plt.subplot2grid((10,4), (0,0), colspan=4,rowspan=8)
	ax2=  plt.subplot2grid((10,4), (8 ,0),colspan=4,rowspan=2,sharex=ax1)
	ax1.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
	ax2.tick_params(axis='y', which='both',left='off',right='on',labelright='on', labelleft='off')
	b_circles=[[matplotlib.patches.Circle((i,3-n),0.20,color='black') for n,c in enumerate(config) if c] for i,config in enumerate([(0,0,1),(0,1,0),(1,0,0)])]
	b_circles=[c for l in b_circles for c in l]
	g_circles=[[matplotlib.patches.Circle((i,3-n),0.20,color='grey') for n,c in enumerate(config) if not c] for i,config in enumerate([(0,0,1),(0,1,0),(1,0,0)])]
	g_circles=[c for l in g_circles for c in l]
	boxes  =[matplotlib.patches.Rectangle((i -0.75/2 ,0), 0.75, 4, lw=4,edgecolor='black',fill=False) for i in range(3)]
	[ax2.add_artist(artist) for artist in boxes]
	[ax2.add_artist(artist) for artist in b_circles]
	[ax2.add_artist(artist) for artist in g_circles]
	ax2.set_yticklabels(["5'",'','N -3','','N -2','','N -1','',"3'"])
	ax2.set_xticklabels(['','N -3','','N -2','','N -1',''])

	ax2.set_ylim([0,4])
	ax2.set_xlim([-0.5,2.5])
	return(fig,(ax1,ax2))




def build_rval_plot(plt):	
	fig = plt.figure()
	ax  = fig.add_subplot(111)
	patches = [matplotlib.patches.Rectangle((0,1), 78, 0.1, facecolor='black'      ),#zrux
	                   matplotlib.patches.Rectangle((281,    1), 22, 0.1, facecolor='red'  ),#uas1
		        matplotlib.patches.Rectangle((389,    1), 22, 0.1, facecolor='red' ), # uass2
	                   matplotlib.patches.Rectangle((551,    1), 6, 0.1, facecolor='green'  ), # TATA
	                   matplotlib.patches.Rectangle((607,    1), 1484, 0.1, facecolor='cyan'  ), # tss 
	                   matplotlib.patches.Rectangle((652,    1), 1402, 0.1, facecolor='blue'  ), # orf
	                   matplotlib.patches.Rectangle((2136,    1), 128, 0.1, facecolor='black'  ), # lexa
	       ]
	[ax.add_patch(p)  for p in patches]
	plt.ylim([0,1.1])
	plt.xlim([0,2246])
	return(fig,ax)


def build_gene_plot(plt):
	fig = plt.figure()
	ax1 = plt.subplot2grid((10,1), (0,0), colspan=1,rowspan=1)
	ax2=  plt.subplot2grid((10,1), (1 ,0),colspan=1,rowspan=9,sharex=ax1)
	ax1.tick_params(axis='y',which='both',labelleft='off',labelright='off',labelbottom='off')
	patches = [matplotlib.patches.Rectangle((48 ,0), 31, 1, facecolor='black'      ),#zrux
	                   matplotlib.patches.Rectangle((281,    0), 21, 1, facecolor='red'  ),#uas1
	                   matplotlib.patches.Rectangle((391,    0), 21, 1, facecolor='red' ), # uass2
	                   matplotlib.patches.Rectangle((551,    0), 6, 1, facecolor='green'  ), # TATA
	                   matplotlib.patches.Rectangle((607,    0), 1484, 1, facecolor='cyan'  ), # tss 
	                   matplotlib.patches.Rectangle((652,    0), 1404, 1, facecolor='blue'  ), # orf
	                   matplotlib.patches.Rectangle((2141,    0), 106, 1, facecolor='black'  ), # lexa
	       ]
	[ax1.add_patch(p)  for p in patches]
	ax1.grid(False)
	ax1.set_ylim([0,1])
	ax1.set_xlim([0,2246])
	return(fig,ax2)



files    =glob.glob("C:\\Users\\Robert\\Desktop\\NewT&ADatabase\\*\\*\\*\\molecule_info.txt")
data     =Analysis(files)
doody    =Subset(data._mols[data._strain=='doody'],color='green',label= 'chd1 tata ACTIVE')
active   =Subset(data._mols[data._strain=='yM8']  ,color='red'  ,label= 'CHD1 TATA ACTIVE')
repressed=Subset(data._mols[data._strain=='yM2']  ,color='blue',label=  'CHD1 TATA repressed')
chd1     =Subset(data._mols[data._strain=='yM208'],color='purple',label='chd1 TATA ACTIVE')
subsets=[active,repressed,chd1,doody]