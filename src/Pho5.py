__version__="01.00.00"
__author__ ="Robert Shelansky"
ZRS =(0,78)
UAS1=(281,303)
UAS2=(389,411)
TATA=(551,557)
TSS = 607
ORF =(652,2055)
NUCLEOSOME_SIZE  =147
NUCLEOSOME_CUTOFF=90

class Configuration:
	"""Configuration is simply a representation of a promoter configuration,
	it holds references to the bubble assigned to the N-1 N-2 and N-3 Positions.
	"""
	#3'-----5'
	#N-1,N-2,N-3
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

	def __init__(self,N_1=None,N_2=None,N_3=None):
		self.n1 = N_1
		self.n2 = N_2
		self.n3 = N_3
		self.configuration = self.IMAP[(    int(bool(self.n1)),
						    				int(bool(self.n2)),
						    				int(bool(self.n3))  )]
		self.size = sum(self.MAP[self.configuration])


	def __str__(self):
		return(str(self.configuration))
  
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


