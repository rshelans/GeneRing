import time
import os.path
import Trace
import scipy
__version__="01.00.00"
__author__ ="Robert Shelansky"
class TraceModel:
	"""
	A model/controller of data which knows how to do specific analysis and maintenance of settings etc.
	Can be used as a data dump for communication with the TraceView.
	"""
	def __init__(self,imgres=scipy.NAN,directory=None,user=None,out_path=None,title=None,version=None,files=None,length=None,threshold=None,smooth=None,**kw):
		self.settings = {
			'length':length,
			'threshold':threshold,
			'smooth':smooth,
			'end':None,
			'midpoint':None}
		self.context  = {
			'files'  :files,
			'version':version,
			'date':time.strftime("%d/%m/%Y %H:%M:%S"),
			'path':None,
			'trace':None,
			'coordinates':None,
			'index':None,
			'domain':None,
			'basename':None,
			'title':title,
			'scale':None,
			'user':user,
			'resolution':None,
			'steadiness':None,
			'saved':['' for i in range(len(files))],
			'out_path':out_path,
			'imgres':imgres,
			'pimgres':None,
			'directory':directory}
		self.molecule = {
			'smoothed':None,
			'fr':None,
			'fl':None,
			'rr':None,
			'rl':None,
			'zipped':None,
			'zlabels':None,
			'symmerty':None,
			'starts':None,
			'bubs':None}

		self.seek(0)
		self.find_midpoint()
		self.analyze()
		

	def find_midpoint(self):
		self.settings['end'        ]=self.context['trace'].edgebuffer(
			threshold=self.settings['threshold'],
			smooth   =self.settings['smooth'])
		self.settings['midpoint'   ]=self.context['trace'].midpoint(
			threshold  =self.settings['threshold'],
			smooth     =self.settings['smooth'],
			end        =self.settings['end'])


	def analyze(self):
		midpoint,smoothed,(fr,fl),(rr,rl),(regions,labels)=self.context['trace'].solve_molecule(
			self.settings['midpoint'],
			self.settings['threshold'],
			self.settings['smooth'],
			self.settings['end'])

		self.molecule['segments'] =len(fr) if len(fr) == len(rr) else float('Nan')
		self.molecule['smoothed'] =smoothed
		self.molecule['fr'      ] =fr
		self.molecule['fl'      ] =fl
		self.molecule['rr'      ] =rr
		self.molecule['rl'      ] =rl
		self.molecule['zipped'  ] =regions
		self.molecule['zlabels' ] =labels
		self.molecule['rlabels' ] =self.context['trace'].label(fl,rl) 
		self.molecule['molecule'] =self.context['trace'].moleculify(
			fr,
			fl,
			rr,
			rl,
			self.settings['length'])
		self.molecule['symmetry'] =self.context['trace'].sd_segment_size         (fr,fl,rr,rl) * self.context['scale']**2
		self.molecule['starts'  ] =self.context['trace'].msd_of_linker_end_points(fr,fl,rr,rl) * self.context['scale']**2
		self.molecule['bubs'    ] =self.context['trace'].msd_of_region_sizes     (fr,fl,rr,rl) * self.context['scale']**2

	def seek(self, index):
		BASE_TO_NM_CONVERSION_FACTOR=0.34#nm/bp
		self.context['path' ]   =self.context['files'][index]
		self.context['index']   =index
		self.context['basename']=os.path.basename(self.context['path'])
		##reads _trace file into  scipy.array
		self.context['coordinates'] = scipy.genfromtxt(
			self.context['path' ],
			delimiter='\t')
		self.context['trace']      = Trace.Trace(self.context['coordinates'])
		self.context['scale']      = self.context['trace'].scale(self.settings['length'])
		self.context['resolution'] = self.context['trace']._ld.mean() * self.context['scale']
		self.context['steadiness'] = scipy.sqrt(self.context['trace']._ld.var()) * self.context['scale']
		self.context['pimgres']    = self.context['scale'] * BASE_TO_NM_CONVERSION_FACTOR
		self.context['domain']     = scipy.array(range(len(self.context['trace'])))

	def write_comments(self,file):
		print("""
##Image_Resolution\t{:.2f} nm/px       
##Predicted_Image_Resolution\t{:>.2f} nm/px
##Tracer\t{}
##Length\t{} bp
##Edgebuffer\t{} Coordinates
##Threshold\t{:.2f} bp
##Smooth\t{} Coordinates
##Midpoint\t{} Coordinate
##Scale\t{:.2f} bp/AU
##Resolution\t{:.2f} bp
##Steadiness\t{:.2f} bp
##Segments\t{:} # 
##Symmetry\t{:.2f} bp
##Linker\t{:.2f} bp
##Region\t{:.2f} bp""".format(
			self.context['imgres'],
			self.context['pimgres'],
			self.context['user'],
			self.settings['length'],
			self.settings['end'],
			self.settings['threshold'] * self.context['scale'],
			self.settings['smooth'],
			self.settings['midpoint'],
			self.context['scale'],
			self.context['resolution'],
			self.context['steadiness'],
			self.molecule['segments'],
			self.molecule['symmetry'],
			self.molecule['starts'],
			self.molecule['bubs']),file=file)

	def save(self):
		base=os.path.basename(self.context['path']).split('.')[0]
		if self.context['out_path'] is not None:
			path=self.context['out_path']
		else:
			path= os.path.dirname(self.context['path'])
			if self.context['directory']:
				path=path+'\\'+base
				if not os.path.exists(path):
					os.makedirs(path)

		mol_file='{}\\{}.mol'.format(path,base)
		reg_file='{}\\{}.reg'.format(path,base)
		with open(mol_file, 'w') as file:
			self.write_comments(file)
			self.molecule['molecule'].write(file)
		with open(reg_file,'w') as file:
			self.write_comments(file)
			reg="\n".join(['{}\t{}\t{}'.format(l,s,e) for (s,e),l in zip(self.molecule['zipped'],self.molecule['zlabels'])])
			print(reg,file=file)
		self.context['saved'][self.context['index']] = 'Saved.'
