import matplotlib
import scipy

__version__="01.00.00"
__author__ ="Robert Shelansky"

ZRS   =(0  , 78)
UAS1  =(281,303)
UAS2  =(389,411)
TATA  =(551,557)
TSS   = 607
ORF   =(652,2055)
LEXA  =(2140,2246)
LENGTH=2246
NUCLEOSOME_SIZE  =147
NUCLEOSOME_CUTOFF=90

class Configuration:
	"""Configuration is simply a representation of a promoter configuration,
	it holds references to the bubble assigned to the N-1 N-2 and N-3 Positions.
	"""
	#		--Pho5--0->
	#		5'--------3'
	#		N-3,N-2,N-1
	MAP={	0:(1,1,1),
			1:(1,1,0),
			2:(1,0,1),
			3:(0,1,1),
			4:(0,0,1),
			5:(0,1,0),
			6:(1,0,0),
			7:(0,0,0)
	}
	IMAP={value:key for key,value in MAP.items()}

	def __init__(self,N_3=None,N_2=None,N_1=None):
		self.n3 = N_3
		self.n2 = N_2
		self.n1 = N_1
		self.configuration = self.IMAP[(    int(bool(self.n3)),
						    				int(bool(self.n2)),
						    				int(bool(self.n1))  )]
		self.size = sum(self.MAP[self.configuration])


	def __str__(self):
		return(str(self.configuration))

	def __len__(self):
		return(sum(self.MAP[self.configuration]))

	def __int__(self):
		return(self.configuration)

	def __eq__(self, other):
		if isinstance(other, Configuration):
			return self.configuration == other.configuration
		elif isinstance(other, (int, float)):
			return self.configuration == other
		else:
			return NotImplemented


#Robert's Regions
N1REGION=[UAS2[1]+24,TSS+NUCLEOSOME_SIZE//2+5]
N3REGION=[ZRS[0],UAS2[0]-35]
N2REGION=[UAS2[0]-83,UAS2[1]+NUCLEOSOME_CUTOFF+10]
#Saketh's Regions
sN1REGION=[490,580]
sN2REGION=[UAS2[0]-83,UAS2[1]+NUCLEOSOME_CUTOFF+10]
sN3REGION=[ZRS[0],UAS2[0]-35]

def _caller_(molecule):
	"""
	 takes a molecule object and returns its estimation of the 
	configuration of the molecuels promoter as defined by 
	
	'Nucleosomal promoter variation generates gene expression noise'
	Christopher R. Brown and Hinrich Boeger1

	The return type is the configuration Number.
	"""

	#Definitions for molecule Regions
	n1   =    [bub for bub in molecule.getOverlap(N1REGION[0],N1REGION[1],NUCLEOSOME_CUTOFF) if bub.isbubble]
	n2   =    [bub for bub in molecule.getExc(UAS2[0],UAS2[1])  if (bub.isbubble and bub.size >NUCLEOSOME_CUTOFF) ]
	n3   =    [bub for bub in molecule.getOverlap(N3REGION[0],N3REGION[1],NUCLEOSOME_CUTOFF) if bub.isbubble]
	return(Configuration(n3,n2,n1))

def bincount_configuration(molecules, config_caller = _caller_):
	return(scipy.bincount(
		scipy.array(list(map(config_caller, molecules))).astype(int),
		minlength= len(Configuration.MAP)))




def build_gene_plot(plt):
	fig    = plt.figure(facecolor='0.8')
	ax2    = plt.axes([0.05,0.1  ,0.9 ,0.85]) 
	ax1    = plt.axes([0.05,0.075,0.9 ,0.025])
	ax1.set_axis_bgcolor('white')
	ax2.set_axis_bgcolor('white')

	ax1.tick_params(axis='y',which='both',top='off',bottom='off',labelleft='off',labelright='off',labelbottom='off')
	ax1.tick_params(axis='x',which='both',top='off',bottom='off')

	ax2.tick_params(axis='x',which='both',labelbottom='off')
	ax2.grid(color='grey',linestyle='--',alpha=0.5)
	ax1.set_xticks(range(-500,2250,250))
	ax2.set_xticks(range(-500,2250,250))
	[i.set_visible(True) for i in ax2.spines.values()]
	[i.set_color('black') for i in ax2.spines.values()]
	[i.set_linewidth(1) for i in   ax2.spines.values()]
	[i.set_visible(True) for i in  ax1.spines.values()]
	[i.set_color('black') for i in ax1.spines.values()]
	[i.set_linewidth(1) for i in   ax1.spines.values()]
	ax2.get_xgridlines()[2].set_linestyle('-')

	patches = [
				matplotlib.patches.Rectangle((0   -TSS ,    0), 30  , 1, facecolor='black'  ), #lexa
				matplotlib.patches.Rectangle((48  -TSS ,    0), 30  , 1, facecolor='dimgrey'), #zrux
				matplotlib.patches.Rectangle((281 -TSS ,    0), 21  , 1, facecolor='red'    ), #uas1
				matplotlib.patches.Rectangle((391 -TSS ,    0), 21  , 1, facecolor='red'    ), #uass2
				matplotlib.patches.Rectangle((551 -TSS ,    0), 6   , 1, facecolor='green'  ), #TATA
				matplotlib.patches.Rectangle((607 -TSS ,    0), 1484, 1, facecolor='cyan'   ), #tss 
				matplotlib.patches.Rectangle((652 -TSS ,    0), 1404, 1, facecolor='blue'   ), #orf
				matplotlib.patches.Rectangle((2140-TSS ,    0), 106 , 1, facecolor='black'  ), #lexa
	]
	[ax1.add_patch(p)  for p in patches]
	#ax1.grid(False)
	ax2.set_ylim([0,1])
	ax2.set_xlim([0-607,2246-607])
	ax1.set_xlim([0-607,2246-607])
	##arrow
	#ax2.plot((0,0)  ,(0, 0.075)    ,linewidth=4, color='black')
	#ax2.plot((0,120),(0.075,0.075)  ,linewidth=4, color='black')
	#ax2.arrow(80,0.075,0.01,0,  color='black',head_starts_at_zero=True,head_length=35)
	return(fig,ax2)

def build_config_plot(plt):
	fig = plt.figure()
	ax1 = plt.subplot2grid((10,4), (0,0), colspan=4,rowspan=8)
	ax2 = plt.subplot2grid((10,4), (8 ,0),colspan=4,rowspan=2,sharex=ax1)
	ax1.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
	ax1.tick_params(axis='y', which='both',right='off')
	ax2.tick_params(axis='y', which='both',left='off',right='off',labelright='on', labelleft='off')
	ax2.tick_params(axis='x',which='both',top='off',bottom='off')
	circles=[[matplotlib.patches.Circle((i,n+1),0.20,color='black') for n,c in enumerate(config) if c] for i,config in Configuration.MAP.items()]
	circles=[c for l in circles for c in l]
	boxes  =[matplotlib.patches.Rectangle((i -0.75/2 ,0), 0.75, 4, lw=4,edgecolor='black',fill=False) for i in range(8)]
	[ax2.add_artist(artist) for artist in boxes]
	[ax2.add_artist(artist) for artist in circles]
	ax2.set_yticklabels(["5'",'','N -3','','N -2','','N -1','',"3'"])
	ax2.set_ylim([0,4])
	ax2.set_xlim([-0.5,7.5])
	return(fig,ax1)

	
def build_nuc_config_plot(plt):	
	fig = plt.figure()
	ax1 = plt.subplot2grid((10,4), (0,0), colspan=4,rowspan=8)
	ax2 = plt.subplot2grid((10,4), (8 ,0),colspan=4,rowspan=2,sharex=ax1)
	ax1.tick_params(axis='y', which='both',right='off')
	ax1.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
	ax2.tick_params(axis='y', which='both',left='off',right='off',labelright='on', labelleft='off')
	ax2.tick_params(axis='x',which='both',top='off',bottom='off')
	b_circles=[[matplotlib.patches.Circle((i,n+1),0.20,color='black') for n,c in enumerate(config) if c] for i,config in enumerate([(1,0,0),(0,1,0),(0,0,1)])]
	b_circles=[c for l in b_circles for c in l]
	g_circles=[[matplotlib.patches.Circle((i,n+1),0.20,color='grey') for n,c in enumerate(config) if not c] for i,config in enumerate([(1,0,0),(0,1,0),(0,0,1)])]
	g_circles=[c for l in g_circles for c in l]
	boxes  =[matplotlib.patches.Rectangle((i -0.75/2 ,0), 0.75, 4, lw=4,edgecolor='black',fill=False) for i in range(3)]
	[ax2.add_artist(artist) for artist in boxes]
	[ax2.add_artist(artist) for artist in b_circles]
	[ax2.add_artist(artist) for artist in g_circles]
	ax2.set_yticklabels(["5'",'','N -3','','N -2','','N -1','',"3'"],weight='bold')
	ax2.set_xticklabels(['','N -3','','N -2','','N -1',''],weight='bold')
	ax2.set_ylim([0,4])
	ax2.set_xlim([-0.5,2.5])
	return(fig,ax1)

