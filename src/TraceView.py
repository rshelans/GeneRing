import matplotlib.widgets
import scipy
__version__="01.00.00"
__author__ ="Robert Shelansky"

def scatter(trace,plt,axes,regions,labels,cmap='bwr'):
	"""
	plots a bubbule and linker coordinates color coded by regions,labels
	"""
	coolabs=scipy.concatenate([[l]*(e-s) for (s,e),l in zip(regions,labels)])
	return(axes.scatter(trace[:,0], trace[:,1],c=coolabs,cmap=plt.get_cmap(cmap) ,lw=0))

def label(labels):
	scipy.random.seed(0)
	rcolors=scipy.array(list(range(max(labels)+1)))
	rcolors=rcolors%2*max(labels)//2+rcolors//2
	rcolors=rcolors/max(rcolors)
	scipy.random.shuffle(rcolors)
	rcolors    =  {i:r for i,r in enumerate(rcolors)}
	cols       =scipy.array([rcolors[i] for i in labels])
	return(labels,cols)

def text_regions(plt,ax,x,y,regions,labels,rlabels,rcols,**kw):
		rcols = plt.get_cmap('rainbow')(rcols)
		t = ax.transData
		t2= ax.transData
		for (s,e),l,r,c in reversed(list(zip(regions,labels,rlabels,rcols))):
			bubcolor=plt.get_cmap('bwr')(float(l))
			region_text = ax.text(
				x,
				y,
				'{:<2d}'.format(r),
				transform=t,
				color=c,
				va='bottom',
				fontdict={'family':'monospace'},
				**kw)
			region_text.draw(ax.figure.canvas.get_renderer())
			t2= matplotlib.transforms.offset_copy(
				region_text._transform,
				x=region_text.get_window_extent().width+10,
				units='dots')
			bubble_text = ax.text(
				x,
				y,
				str(l),
				transform=t2,
				color=bubcolor,
				va='bottom',
				fontdict={'family':'monospace'},
				**kw)	
			t2= matplotlib.transforms.offset_copy(
				region_text._transform,
				x=region_text.get_window_extent().width+20,
				units='dots')
			string= '{:>6d}..{:>5d}  {:>4d} bp'.format(s,e,e-s)
			coord_text = ax.text(
				x,
				y,
				string,
				transform=t2,
				color='black',
				va='bottom',
				fontdict={'family':'monospace'},
				**kw)
			t=matplotlib.transforms.offset_copy(
				region_text._transform,
				y=region_text.get_window_extent().height,
				units='dots')

def label_scatter (axes, coordinates, labels, colors):
		"""
		labels each region with some label in a color 
		"""
		txts=[]
		for (x,y),l,c in zip(coordinates,labels,colors):
			txt=axes.text(x,y, l,color='black',bbox=dict(facecolor=c))
			txts.append(txt)
		return(txts)
class TraceView:
	"""
	builds a window which can show a trace and adjust parameter settings .
	Builds a GUI.
	"""
	def __init__(self, plt, model):
		self.plt     = plt
		self.model   = model
		self.fig     = self.plt.figure(facecolor='0.8')
		#Panel for Distance to closest point
		self.panel1  = self.plt.axes([0.05,0.15 ,0.7 ,0.1])  
		#Panels for trace
		self.panel2  = self.plt.axes([0.05,0.30 ,0.35,0.6])
		self.panel3  = self.plt.axes([0.40,0.30 ,0.35,0.6])
		#textPanel    
		self.panel4  = self.plt.axes([0.755,0.3,0.24,0.7],axisbg='0.8')
		self.panel4.set_xlim([0,1])
		self.panel4.set_ylim([0,1])
		self.panel4.axis('off')
		#Slider Panels
		self.wdj1    = self.plt.axes([0.10,0.025,0.60,0.02])
		self.wdj2    = self.plt.axes([0.10,0.050,0.60,0.02])
		self.wdj3    = self.plt.axes([0.10,0.075,0.60,0.02])
		self.wdj8    = self.plt.axes([0.10,0.1  ,0.60,0.02])
		#Panel for buttons
		self.wdj4    = self.plt.axes([0.8 ,0.025,0.15,0.02])
		self.wdj5    = self.plt.axes([0.8 ,0.050,0.15,0.02])
		self.wdj6    = self.plt.axes([0.8 ,0.075,0.15,0.02])
		self.wdj7    = self.plt.axes([0.8 ,0.1  ,0.15,0.02])

		self.saved   = self.fig.text(0.875, 0.02,'', va='top',ha='center') 

		self.mp_slider = matplotlib.widgets.Slider(
			self.wdj1,
			'midpoint',
			1,
			len(self.model.context['trace']),
			valinit =self.model.settings['midpoint'],
			dragging=False,
			valfmt  ='%d')
		self.td_slider = matplotlib.widgets.Slider(
			self.wdj2,
			'threshold',
			0.1,
			50,
			valinit =self.model.settings['threshold'],
			dragging=False,
			valfmt  ='%.2f')
		self.sy_slider = matplotlib.widgets.Slider(
			self.wdj3,
			'smooth',
			1,
			50,
			valinit = self.model.settings['smooth'],
			dragging=False,
			valfmt  ='%d')

		self.se_slider=matplotlib.widgets.Slider(
			self.wdj8,
			'seek',
			0,
			len(self.model.context['files']),
			valinit = self.model.context['index'],
			dragging=False,
			valfmt  ='%d'
			)

		self.sv_button = matplotlib.widgets.Button(self.wdj4,'SAVE')
		self.nx_button = matplotlib.widgets.Button(self.wdj5,'>> NEXT >>')
		self.pv_button = matplotlib.widgets.Button(self.wdj6,'<< PREV <<')
		self.mp_button = matplotlib.widgets.Button(self.wdj7,'Find Midpoint')


		self.panel1.set_xlim([0, len(self.model.context['trace'])-1])
		self.panel2.tick_params(
						axis='both',
						direction='out',
						which='both',
						bottom='off',
						labelbottom='off',
						left='off',
						labelleft='off',
						right='off',
						labelright='off',
						top='off',
						labeltop='off')
		self.panel3.tick_params(
						axis='both',
						direction='out',
						which='both',
						bottom='off',
						labelbottom='off',
						left='off',
						labelleft='off',
						right='off',
						labelright='off',
						top='off',
						labeltop='off')

		self.mid_plot       = self.panel1.axvline(
			x=self.model.settings['midpoint']-1,
			color='red' )
		self.thresh_plot    = self.panel1.axhline(
			y=self.model.settings['threshold'] ,
			color='cyan')
		self.line_plot,     = self.panel1.plot(
			self.model.context['domain'], 
			self.model.molecule['smoothed'])	
		self.mol_plot       = scatter(
			self.model.context['trace'],
			self.plt,
			self.panel2,
			self.model.molecule['zipped'],
			self.model.molecule['zlabels'])
		self.mpf_plot       = self.panel2.scatter(
			self.model.context['trace'][self.model.settings['midpoint'],0],
			self.model.context['trace'][self.model.settings['midpoint'],1],
			color='black',
			s=100)
		self.mpr_plot       = self.panel2.scatter(
			self.model.context['trace'][self.model.settings['midpoint']+1,0],
			self.model.context['trace'][self.model.settings['midpoint']+1,1],
			color='black',
			s=100)

		reg_labs, reg_cols =label(self.model.molecule['rlabels'])
		self.reg_plot      =scatter(
			self.model.context['trace'],
			self.plt,
			self.panel3,
			self.model.molecule['zipped'],
			reg_cols,
			cmap='rainbow')

		regs  = self.model.molecule['fr'] if len(self.model.molecule['fr']) > len(self.model.molecule['rr']) else self.model.molecule['rr']
		mids  = [self.model.context['trace']._midpoint(s,e-1) for s,e in regs]
		coords= [self.model.context['trace'][mid] for mid in mids]
		labs  = [i            for i   in range(len(regs))]
		_,cols= self.plt.get_cmap('rainbow')(label(labs))
		self.lab_txts =label_scatter(
			self.panel3,
			coords,
			labs,
			cols)
		
		self.title =self.fig.text(
			0.05,
			1,
			self.build_title(),
			va='top',
			ha='left',
			fontdict={'family':'monospace'})
		report     =self.build_report()
		self.report=self.panel4.text(
			0,
			1,
			report,
			va='top',
			fontdict={'family':'monospace'},
			fontweight='bold')
		####
		##REGION TEXT
		####
		
		mol=self.model.molecule['molecule']
		if mol is not None:
			regs=[[reg.start,reg.end] for reg in mol     ]
			bubs=[reg.isbubble for reg in mol            ]
			labs=list(reversed([i            for i   in range(len(regs))]))
			_,cols=label(labs)
			self.regtext=text_regions(
				self.plt,
				self.panel4,
				0,
				0,
				regs,
				bubs,
				labs,
				cols,
				fontweight='bold')
		##REGIONTEXT

		self.mp_slider.on_changed(self.toggle_settings)
		self.td_slider.on_changed(self.toggle_settings)
		self.sy_slider.on_changed(self.toggle_settings)
		self.se_slider.on_changed(self.seek  )
		self.sv_button.on_clicked(self.save  )
		self.nx_button.on_clicked(self.next  )
		self.pv_button.on_clicked(self.prev  )
		self.mp_button.on_clicked(self.find_midpoint)




	def build_title(self):
		FILE_INFO="""
{}
Version : v{}
Path    : {}
File    : {}
Date    : {}
""".format(
		self.model.context['title'],
		self.model.context['version'],
		self.model.context['path'    ],
		self.model.context['basename'],
		self.model.context['date'])
		return(FILE_INFO)

	def build_report(self):
		SETTINGS_INFO="""
SETTINGS
Length    : {} bp
EdgeIgnore: {} Coordinates
Threshold : {:0.2f} bp
Smooth    : {} Coordinates
Midpoint  : {} Coordinate
""".format(
		self.model.settings['length'],
		self.model.settings['end'],
		self.model.settings['threshold'] * self.model.context['scale'],
		self.model.settings['smooth'],
		self.model.settings['midpoint'])

		TRACE_INFO="""
TRACE INFO
Image Resolution: {:>10.2f} nm/px       
Pred. Image Res.: {:>10.2f} nm/px
Trace Scale     : {:>10.2f} bp/AU
Trace Resolution: {:>10.2f} bp
Trace Steadiness: {:>10.2f} bp
""".format(
		self.model.context['imgres'],
		self.model.context['pimgres'],
		self.model.context['scale'],
		self.model.context['resolution'],
		self.model.context['steadiness'])

		if self.model.molecule['symmetry'] is not None:
			TRACE_FIT="""
BUBBLE CALL INFO
Segments      : {:>10} # 
Symmetry score: {:>10.2f} bp²
Linker score  : {:>10.2f} bp²
Region score  : {:>10.2f} bp²
""".format(
		self.model.molecule['segments'],
		self.model.molecule['symmetry'],
		self.model.molecule['starts'],
		self.model.molecule['bubs'])
		else:
			TRACE_FIT="""
COULD NOT ANALYZE TRACE
    --try alternate settings
	"""
		return(SETTINGS_INFO+TRACE_INFO+TRACE_FIT)


	def update(self):	
		self.saved.set_text(self.model.context['saved'][self.model.context['index']])
		self.panel1.set_ylim  ([0,max(self.model.molecule['smoothed'])])
		self.panel1.set_xlim  ([0,len(self.model.molecule['smoothed'])])
		self.line_plot.set_ydata(self.model.molecule['smoothed'] )
		self.line_plot.set_xdata(self.model.context['domain'])
		self.mid_plot.set_xdata    (self.model.settings['midpoint'] )
		self.thresh_plot.set_ydata (self.model.settings['threshold'])
		
		reg_labs, reg_cols=label(self.model.molecule['rlabels'])
		reg_cols =[[l]*(e-s) for (s,e),l in zip(self.model.molecule['zipped'],reg_cols)]
		reg_cols =self.plt.get_cmap('rainbow')(scipy.concatenate(reg_cols).astype(float))

		self.reg_plot.set_color(reg_cols)
		self.reg_plot.set_offsets((self.model.context['trace']))
		offset=40
		self.panel3.set_ylim([min(self.model.context['trace'][:,1])-offset,max(self.model.context['trace'][:,1])+offset] )
		self.panel3.set_xlim([min(self.model.context['trace'][:,0])-offset,max(self.model.context['trace'][:,0])+offset] )

		cmap=self.plt.cm.get_cmap('bwr')
		coo_cols=[[l]*(e-s) for (s,e),l in zip(self.model.molecule['zipped'],self.model.molecule['zlabels'])]
		coo_cols=cmap(scipy.concatenate(coo_cols).astype(float))
		self.mol_plot.set_color    (coo_cols)
		self.mol_plot.set_offsets( self.model.context['trace'])

		self.mpf_plot.set_offsets( (self.model.context['trace'] [self.model.settings['midpoint']  ,0],
			self.model.context['trace'][self.model.settings['midpoint']  ,1])) 
		self.mpr_plot.set_offsets( (self.model.context['trace'] [self.model.settings['midpoint']+1,0],
			self.model.context['trace'][self.model.settings['midpoint']+1,1]))
		self.panel2.set_ylim([min(self.model.context['trace'][:,1])-offset,max(self.model.context['trace'][:,1])+offset])
		self.panel2.set_xlim([min(self.model.context['trace'][:,0])-offset,max(self.model.context['trace'][:,0])+offset])

		self.title.set_text(self.build_title())

		self.panel4.cla()
		self.panel4.set_xlim([0,1])
		self.panel4.set_ylim([0,1])
		self.panel4.axis('off')
		report     =self.build_report()
		self.report=self.panel4.text(
			0,
			1,
			report,
			va='top',
			fontdict={'family':'monospace'})

		[txt.remove() for txt in self.lab_txts]
		regs  = self.model.molecule['fr'] if len(self.model.molecule['fr']) > len(self.model.molecule['rr']) else self.model.molecule['rr']
		mids  = [self.model.context['trace']._midpoint(s,e-1) for s,e in regs]
		coords= [self.model.context['trace'][mid] for mid in mids]
		labs  = [i            for i   in range(len(regs))]
		_,cols= self.plt.get_cmap('rainbow')(label(labs))
		self.lab_txts =label_scatter(
			self.panel3,
			coords,
			labs,
			cols)
		##############################
		############Region txt########
		##############################
		mol=self.model.molecule['molecule']
		if mol is not None:
			regs=[[reg.start,reg.end] for reg in mol     ]
			bubs=[reg.isbubble for reg in mol            ]
			labs=list(reversed([i            for i   in range(len(regs))]))
			_,cols=label(labs)
			self.regtext=text_regions(
				self.plt,
				self.panel4,
				0,
				0,
				regs,
				bubs,
				labs,
				cols,
				fontweight='bold')
		###################################
	
	def toggle_settings(self, val):
		self.model.settings['threshold'  ] =float(self.td_slider.val)
		self.model.settings['smooth'] =int  (self.sy_slider.val)
		self.model.settings['midpoint'   ] =int  (self.mp_slider.val) - 1
		self.model.analyze()
		self.update()
		self.fig.canvas.draw_idle()



	def reset_settings(self):
		#reset Model Settings
		self.model.settings['threshold'  ]=self.td_slider.valinit
		self.model.settings['smooth']=self.sy_slider.valinit
		self.model.find_midpoint()
		#reset view settings
		self.mp_slider.val=self.model.settings['midpoint'   ]
		self.td_slider.val=self.model.settings['threshold'  ]
		self.sy_slider.val=self.model.settings['smooth']
		self.mp_slider.set_val(self.model.settings['midpoint'   ])
		self.td_slider.set_val(self.model.settings['threshold'  ])
		self.sy_slider.set_val(self.model.settings['smooth'])

	def seek(self, val):
		self.model.seek(int(self.se_slider.val))
		self.reset_settings()

	def next(self, val):
		next= min(self.model.context['index']+1,len(self.model.context['files'])-1)
		self.se_slider.set_val(next)

	def prev(self,val):
		prev= max(self.model.context['index']-1,0)
		self.se_slider.set_val(prev)

	def save(self,val):
		if self.model.molecule['molecule']:
			self.saved.set_text('Saved.')
			self.model.save()
		else:
			self.saved.set_text('Failed.')

	def find_midpoint(self,val):
		self.model.find_midpoint()
		self.mp_slider.val=self.model.settings['midpoint'   ]
		self.mp_slider.set_val(self.model.settings['midpoint'   ])


	def show(self):
		self.plt.show()
