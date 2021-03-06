import Chromatin
import Molecules
import Pho5

import scipy
import scipy.stats
import matplotlib.pyplot as plt

## This File Contains all kinds of example analysis for looking
## at groups of Molecules using the Chromatin Package. Check it out woot.



PATH_TO_DATABASE = 'E:\\database\\'
data             = Molecules.DataBase(PATH_TO_DATABASE)

s1    = data._mols[data._strain=='yM19']
s2    = data._mols[data._strain=='yM8']
s3    = data._mols[data._strain=='yC169']
s4    = data._mols[data._strain=='yM208']
s5    = data._mols[data._strain=='yM2']
s6    = data._mols[data._strain=='yM89']

# print(data._name[data._check == False])


alpha     =0.95
iterations=1000
x         =list(range(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS))
binwidth  =20






# ##RVALUE ANALYSIS 
# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('R-Value', fontsize=20)
# ax.plot(x, Molecules.rvalue(s1),color='gold'    ,linewidth=2)
# ax.plot(x, Molecules.rvalue(s2),color='darkblue',linewidth=2)

# ##RVALYE ANALYSIS THAT HAS BEEN SMOOTHED
# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('R-Value Smoothed: {}'.format(binwidth), fontsize=20)
# ax.plot(x, Molecules.smooth(Molecules.rvalue(s1),binwidth),color='gold'    ,linewidth=2)
# ax.plot(x, Molecules.smooth(Molecules.rvalue(s2),binwidth),color='darkblue',linewidth=2)

##RVALUE ANALYSIS USING THE BOOTHSTRAP TOOL
# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('R-Value Bootstrapped: {}X'.format(iterations), fontsize=20)

# mew,error,_=Molecules.bootstrap(s1 , Molecules.rvalue, iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='gold',linewidth=0.0)
# ax.plot        (x, mew,color='gold',linewidth=2)

# mew,error,_=Molecules.bootstrap(s2 , Molecules.rvalue, iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='darkblue',linewidth=0.0)
# ax.plot        (x, mew,color='darkblue',linewidth=2)

# mew,error,_=Molecules.bootstrap(s3 , Molecules.rvalue, iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='green',linewidth=0.0)
# ax.plot        (x, mew,color='green',linewidth=2)

# mew,error,_=Molecules.bootstrap(s4 , Molecules.rvalue, iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='purple',linewidth=0.0)
# ax.plot        (x, mew,color='purple',linewidth=2)

# mew,error,_=Molecules.bootstrap(s5 , Molecules.rvalue, iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='red',linewidth=0.0)
# ax.plot        (x, mew,color='red',linewidth=2)

# mew,error,_=Molecules.bootstrap(s6 , Molecules.rvalue, iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='orange',linewidth=0.0)
# ax.plot        (x, mew,color='orange',linewidth=2)



# ##RVALUE SMOOTHED ANALYSIS USING THE BOOT STRAP TOOL
# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('R-Value Bootstrapped & Smoothed: {}X, {}'.format(iterations,binwidth) , fontsize=20)
# mew,error,_=Molecules.bootstrap(s1 , lambda mols: Molecules.smooth(Molecules.rvalue(mols),binwidth), iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='gold',linewidth=0.0)
# ax.plot        (x, mew,color='gold',linewidth=2)

# mew,error,_=Molecules.bootstrap(s2 ,lambda mols: Molecules.smooth(Molecules.rvalue(mols),binwidth), iterations)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='darkblue',linewidth=0.0)
# ax.plot        (x, mew,color='darkblue',linewidth=2)

# #MIDPOINT ANALYSIS
# ##THIS is more complicated though you could simply use the midpoint method on the mols lists wihtout all the fanciness
# ##THIS PARTICULAR ONE LOOKS AT LINKER MIDPOINTS AND then take 1-the smoothed value to get NUC Positions
# binwidth=147/2.5
# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('Midpoint', fontsize=20)



# mew,error,_=Molecules.bootstrap(s1,lambda mols: Molecules.smooth(Molecules.midpoint(mols, gettr=Chromatin.Molecule.getLinkers),binwidth), iterations)  
# ax.fill_between(x, 1-error[0,:], 1-error[1,:], alpha=0.30,color='gold',linewidth=0.0)
# ax.plot   (x, 1-mew,color='gold',linewidth=2)

# mew,error,_=Molecules.bootstrap(s2,lambda mols: Molecules.smooth(Molecules.midpoint(mols, gettr=Chromatin.Molecule.getLinkers),binwidth), iterations)  
# ax.fill_between(x, 1-error[0,:], 1-error[1,:], alpha=0.30,color='darkblue',linewidth=0.0)
# ax.plot   (x, 1-mew,color='darkblue',linewidth=2)

# mew,error,_=Molecules.bootstrap(s5,lambda mols: Molecules.smooth(Molecules.midpoint(mols,gettr=Chromatin.Molecule.getLinkers),binwidth), iterations)  
# ax.fill_between(x, 1-error[0,:], 1-error[1,:], alpha=0.20,color='red',linewidth=0.0)
# ax.plot   (x, 1-mew,color='red',linewidth=2)

# mew,error,_=Molecules.bootstrap(s3,lambda mols: Molecules.smooth(Molecules.midpoint(mols, gettr=Chromatin.Molecule.getLinkers),binwidth), iterations)  
# ax.fill_between(x, 1-error[0,:], 1-error[1,:], alpha=0.30,color='green',linewidth=0.0)
# ax.plot   (x, 1-mew,color='green',linewidth=2)

# mew,error,_=Molecules.bootstrap(s4,lambda mols: Molecules.smooth(Molecules.midpoint(mols,gettr=Chromatin.Molecule.getLinkers),binwidth), iterations)  
# ax.fill_between(x, 1-error[0,:], 1-error[1,:], alpha=0.20,color='purple',linewidth=0.0)
# ax.plot   (x, 1-mew,color='purple',linewidth=2)

# mew,error,_=Molecules.bootstrap(s6,lambda mols: Molecules.smooth(Molecules.midpoint(mols,gettr=Chromatin.Molecule.getLinkers),binwidth), iterations)  
# ax.fill_between(x, 1-error[0,:], 1-error[1,:], alpha=0.20,color='orange',linewidth=0.0)
# ax.plot   (x, 1-mew,color='orange',linewidth=2)


##CONFIGURATION ANALYSIS CONFIGURATIONS
##List of configurations  scipy.array(list(map(Pho5.get_configuration,s1))).astype(int)
##Number of Possible configurations
##List where index is the configuration and the value is the number of tiems that configurations appears
# width=0.35
# x    =scipy.array(range(len(Pho5.Configuration.MAP)))
# s1_configs = Pho5.bincount_configuration(s1)
# s2_configs = Pho5.bincount_configuration(s2)

# fig, ax = Pho5.build_config_plot(plt)
# fig.suptitle('Configuration', fontsize=20)
# ax.bar(x-width, s1_configs/sum(s1_configs), width, color='gold')#, yerr=men_std)
# ax.bar(x      , s2_configs/sum(s2_configs), width, color='darkblue')

##CONFIGURATION ANALYSIS NUCLEOSOME NUMBER
#width=0.35
#x    = scipy.array(range(len([n1,n2,n3])))
#n1   = [0,2,3,4]
#n2   = [0,1,3,5]
#n3   = [0,1,2,6]
#s1_nucs= scipy.array([sum(s1_configs[n1]), sum(s1_configs[n2]), sum(s1_configs[n3])])
#s2_nucs= scipy.array([sum(s2_configs[n1]), sum(s2_configs[n2]), sum(s2_configs[n3])]) 

# fig,ax = Pho5.build_nuc_config_plot(plt)
# fig.suptitle('Configuration Nucleosome Counts', fontsize=20)
# ax.bar(x-width, s1_nucs/len(s1), width, color='gold')#, yerr=men_std)
# ax.bar(x      , s2_nucs/len(s2), width, color='darkblue')

##Bootstrapped Configuration Analysis
# width=0.1
# x    =scipy.array(range(len(Pho5.Configuration.MAP)))
# fig,ax = Pho5.build_config_plot(plt)
# fig.suptitle('Bootstrapped Configuration Analysis', fontsize=20)
# error_config={'ecolor':'black','capsize':0,'elinewidth':2}

# mew,err,_=Molecules.bootstrap(s1, Pho5.bincount_configuration, iterations, alpha)
# ax.bar(x-2*width, mew/len(s1), width, color='gold', yerr=scipy.absolute(mew-err)/len(s1),error_kw=error_config)

# mew,err,_=Molecules.bootstrap(s2, Pho5.bincount_configuration, iterations, alpha)
# ax.bar(x-1*width  , mew/len(s2), width, color='darkblue', yerr=scipy.absolute(mew-err)/len(s2),error_kw=error_config)

# mew,err,_=Molecules.bootstrap(s3, Pho5.bincount_configuration, iterations, alpha)
# ax.bar(x          , mew/len(s3), width, color='green', yerr=scipy.absolute(mew-err)/len(s3),error_kw=error_config)

# mew,err,_=Molecules.bootstrap(s4, Pho5.bincount_configuration, iterations, alpha)
# ax.bar(x+1*width     , mew/len(s4), width, color='purple', yerr=scipy.absolute(mew-err)/len(s4),error_kw=error_config)

# mew,err,_=Molecules.bootstrap(s5, Pho5.bincount_configuration, iterations, alpha)
# ax.bar(x+1*width   , mew/len(s5), width, color='red', yerr=scipy.absolute(mew-err)/len(s5), error_kw=error_config)

# mew,err,_=Molecules.bootstrap(s6, Pho5.bincount_configuration, iterations, alpha)
# ax.bar(x+2*width   , mew/len(s6), width, color='orange', yerr=scipy.absolute(mew-err)/len(s6), error_kw=error_config)





##BOOTSTRAPPED CONFIGURATION NUCLEOSOME NUMBER
# fig,ax   = Pho5.build_nuc_config_plot(plt)
# fig.suptitle('Configuration Nucleosome Counts', fontsize=20)
# width=0.1
# n1   = [0,2,3,4]
# n2   = [0,1,3,5]
# n3   = [0,1,2,6]
# x    = scipy.array(range(len([n1,n2,n3])))
# error_config={'ecolor':'black','capsize':0,'elinewidth':2}

# def GETNUCS(molecules):
# 	n1   = [0,2,3,4]
# 	n2   = [0,1,3,5]
# 	n3   = [0,1,2,6]
# 	confs =Pho5.bincount_configuration(molecules)
# 	return(scipy.array([sum(confs[n3]), sum(confs[n2]), sum(confs[n1])]))

# mew,err,_=Molecules.bootstrap(s1, GETNUCS, iterations, alpha)
# ax.bar(x-3*width, mew/len(s1), width, color='gold',yerr=scipy.absolute(mew-err)/len(s1),error_kw=error_config)#, yerr=men_std)

# mew,err,_=Molecules.bootstrap(s6, GETNUCS, iterations, alpha)
# ax.bar(x-2*width, mew/len(s6), width, color='orange',yerr=scipy.absolute(mew-err)/len(s6),error_kw=error_config)#, yerr=men_std)

# mew,err,_=Molecules.bootstrap(s2, GETNUCS, iterations, alpha)
# ax.bar(x-1*width,mew/len(s2), width, color='darkblue',yerr=scipy.absolute(mew-err)/len(s2),error_kw=error_config)#, yerr=men_std)

# mew,err,_=Molecules.bootstrap(s3, GETNUCS, iterations, alpha)
# ax.bar(x-0*width, mew/len(s3), width, color='green',yerr=scipy.absolute(mew-err)/len(s3),error_kw=error_config)#, yerr=men_std)

# mew,err,_=Molecules.bootstrap(s4, GETNUCS, iterations, alpha)
# ax.bar(x+1*width, mew/len(s4), width, color='purple',yerr=scipy.absolute(mew-err)/len(s4),error_kw=error_config)#, yerr=men_std)

# mew,err,_=Molecules.bootstrap(s5, GETNUCS, iterations, alpha)
# ax.bar(x+2*width, mew/len(s5), width, color='red',yerr=scipy.absolute(mew-err)/len(s5),error_kw=error_config)#, yerr=men_std)



##HEATMAP OF NUCLEOSOMES
# EACH NUCLEOSOME VS Every other Nuclesome
# fig = plt.figure()
# fig.suptitle('Configuration Nucleosome Counts <TATA>', fontsize=20)
# ax = plt.subplot(111)
# confs = list( map(Pho5._caller_, s1) ) 
# confs = scipy.array([Pho5.Configuration.MAP[int(c)] for c in confs])
# ax.imshow(scipy.corrcoef(confs.transpose()),interpolation='none',vmin=0,vmax=1,cmap='cool')

# fig = plt.figure()
# fig.suptitle('Configuration Nucleosome Counts <TATA>(Pho4D85-99)', fontsize=20)
# ax = plt.subplot(111)
# confs = list( map(Pho5._caller_, s6) ) 
# confs = scipy.array([Pho5.Configuration.MAP[int(c)] for c in confs])
# ax.imshow(scipy.corrcoef(confs.transpose()),interpolation='none',vmin=0,vmax=1,cmap='cool')


# fig = plt.figure()
# fig.suptitle('Configuration Nucleosome Counts ACTIVE', fontsize=20)
# ax = plt.subplot(111)
# confs = list( map(Pho5._caller_, s2) ) 
# confs = scipy.array([Pho5.Configuration.MAP[int(c)] for c in confs])
# ax.imshow(scipy.corrcoef(confs.transpose()),interpolation='none',vmin=0,vmax=1,cmap='cool')

# fig = plt.figure()
# fig.suptitle('Configuration Nucleosome Counts <chd1><tata>', fontsize=20)
# ax = plt.subplot(111)
# confs = list( map(Pho5._caller_, s3) ) 
# confs = scipy.array([Pho5.Configuration.MAP[int(c)] for c in confs])
# ax.imshow(scipy.corrcoef(confs.transpose()),interpolation='none',vmin=0,vmax=1,cmap='cool')

# fig = plt.figure()
# fig.suptitle('Configuration Nucleosome Counts <cbd1>', fontsize=20)
# ax = plt.subplot(111)
# confs = list( map(Pho5._caller_, s4) ) 
# confs = scipy.array([Pho5.Configuration.MAP[int(c)] for c in confs])
# ax.imshow(scipy.corrcoef(confs.transpose()),interpolation='none',vmin=0,vmax=1,cmap='cool')

# fig = plt.figure()
# fig.suptitle('Configuration Nucleosome Counts REPRESSED', fontsize=20)
# ax = plt.subplot(111)
# confs = list( map(Pho5._caller_, s5) ) 
# confs = scipy.array([Pho5.Configuration.MAP[int(c)] for c in confs])
# ax.imshow(scipy.corrcoef(confs.transpose()),interpolation='none',vmin=0,vmax=1,cmap='cool')


# #HEARMAP OF MOLECULE
# fig, ax = Pho5.build_gene_plot(plt)
# extent  = [-Pho5.TSS,Pho5.LENGTH-Pho5.TSS,-Pho5.TSS,Pho5.LENGTH-Pho5.TSS]
# fig.suptitle('Full Molecule HEATMAP', fontsize=20)
# c=scipy.corrcoef(scipy.array([mol._array_ for mol in s6]).transpose())
# ax.imshow(c, extent=extent, origin='lower', cmap='bwr', vmin=-1, vmax=1,interpolation='none')
# ax.set_ylim([-Pho5.TSS,Pho5.LENGTH-Pho5.TSS])
# ax.grid(False)

##BUBBLE SIZES DISTRIBUTIONS
# fig, ax = Pho5.build_gene_plot(plt)
# fig.suptitle('Bubble Size Distributions', fontsize=20)
# bins  =scipy.append(scipy.arange(0,2225,125),Pho5.LENGTH)
# positions=scipy.array([s+(e-s)/2-Pho5.TSS for s,e in zip(bins[:-1],bins[1:])])

# groups=[list(Molecules.get_bubbles(s1,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True,positions=positions-15, widths=20, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)

# edit  = plt.setp(boxp['boxes']   , edgecolor='gold',facecolor='gold',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='gold'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

# groups=[list(Molecules.get_bubbles(s2,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True, positions=positions+15, widths=20, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)
# ax.set_xticks(positions)
# ax.grid(color='grey',linestyle='--',alpha=0.5)

# edit  = plt.setp(boxp['boxes']   , edgecolor='darkblue',facecolor='darkblue',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='darkblue'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

##ALL LINKE SIZE SIZES DISTRIBUTIONS STEM LEAF
# widths=10
# gap   =5

# groups=[list(Molecules.get_bubbles(s1,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True,positions=positions-2*widths-2*gap, widths=widths, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)

# edit  = plt.setp(boxp['boxes']   , edgecolor='gold' ,facecolor='gold' ,linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='gold'     ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

# groups=[list(Molecules.get_bubbles(s2,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True, positions=positions-1*widths-1*gap, widths=widths, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)
# ax.set_xticks(positions)
# ax.grid(color='grey',linestyle='--',alpha=0.5)

# edit  = plt.setp(boxp['boxes']   , edgecolor='darkblue',facecolor='darkblue',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='darkblue'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

# groups=[list(Molecules.get_bubbles(s3,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True, positions=positions, widths=widths, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)
# ax.set_xticks(positions)
# ax.grid(color='grey',linestyle='--',alpha=0.5)

# edit  = plt.setp(boxp['boxes']   , edgecolor='green',facecolor='green',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='green'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

# groups=[list(Molecules.get_bubbles(s6,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True, positions=positions+1*widths+1*gap, widths=widths, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)
# ax.set_xticks(positions)
# ax.grid(color='grey',linestyle='--',alpha=0.5)

# edit  = plt.setp(boxp['boxes']   , edgecolor='orange',facecolor='orange',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='orange'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

# groups=[list(Molecules.get_linkers(s4,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True, positions=positions+1*widths+1*gap, widths=widths, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)
# ax.set_xticks(positions)
# ax.grid(color='grey',linestyle='--',alpha=0.5)

# edit  = plt.setp(boxp['boxes']   , edgecolor='purple',facecolor='purple',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='purple'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

# groups=[list(Molecules.get_bubbles(s5,s,e)) for s,e in zip(bins[:-1],bins[1:])]
# sizes =[list(map(len,bubs)) for bubs in groups]
# boxp  = ax.boxplot(sizes, notch=True, positions=positions+2*widths+2*gap, widths=widths, patch_artist=True)
# ax.set_ylim(0,max((s for size in sizes for s in size)))
# ax.set_xlim(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS)
# ax.set_xticks(positions)
# ax.grid(color='grey',linestyle='--',alpha=0.5)

# edit  = plt.setp(boxp['boxes']   , edgecolor='red',facecolor='red',linewidth=2)
# edit  = plt.setp(boxp['whiskers'], color='red'    ,linestyle='-'    ,linewidth=2) 
# edit  = plt.setp(boxp['medians'] , color='white'    ,linestyle='-'    ,linewidth=2)
# edit  = plt.setp(boxp['caps']    , linewidth=0)
# edit  = plt.setp(boxp['fliers']  , markersize=0)

##Bubble SIZE RVALUE STYLE PLOT
# ##CHATT ABOUT MEDIAN AND MEAN MAYBE HERE
# fig, ax = Pho5.build_gene_plot(plt)
# fig.suptitle('Bubble Size By Coordinate Mean', fontsize=20)
# alpha     =0.95
# iterations=1000
# x         =list(range(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS))

# bubsize=scipy.array([list(s) for s in Molecules.linker_size(s1)])
# mew,error,_=Molecules.bootstrap(bubsize,
# 	lambda x: scipy.stats.nanmedian(x,axis=0),
# 	iterations,
# 	alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='gold',linewidth=0.0)
# ax.plot(x,mew,color='gold')

# bubsize=scipy.array([list(s) for s in Molecules.bubble_size(s2)])
# mew,error,_=Molecules.bootstrap(bubsize,
# 	lambda x: scipy.stats.nanmedian(x,axis=0),
# 	iterations,
# 	alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='darkblue',linewidth=0.0)
# ax.plot(x,mew,color='darkblue')

# bubsize=scipy.array([list(s) for s in Molecules.bubble_size(s3)])
# mew,error,_=Molecules.bootstrap(bubsize,
# 	lambda x: scipy.stats.nanmedian(x,axis=0),
# 	iterations,
# 	alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='green',linewidth=0.0)
# ax.plot(x,mew,color='green')

# bubsize=scipy.array([list(s) for s in Molecules.linker_size(s4)])
# mew,error,_=Molecules.bootstrap(bubsize,
# 	lambda x: scipy.stats.nanmedian(x,axis=0),
# 	iterations,
# 	alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='purple',linewidth=0.0)
# ax.plot(x,mew,color='purple')

# bubsize=scipy.array([list(s) for s in Molecules.linker_size(s5)])
# mew,error,_=Molecules.bootstrap(bubsize,
# 	lambda x: scipy.stats.nanmedian(x,axis=0),
# 	iterations,
# 	alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='red',linewidth=0.0)
# ax.plot(x,mew,color='red')

# bubsize=scipy.array([list(s) for s in Molecules.linker_size(s6)])
# mew,error,_=Molecules.bootstrap(bubsize,
# 	lambda x: scipy.stats.nanmedian(x,axis=0),
# 	iterations,
# 	alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.20,color='orange',linewidth=0.0)
# ax.plot(x,mew,color='orange')


# ax.set_ylim([0, scipy.nanmax(bubsize)])



##CORRELATION FROM PROMOTER PLOT
# fig, ax = Pho5.build_gene_plot(plt)
# fig.suptitle('CORRELATION PLOTS', fontsize=20)

# promoter =Pho5.N2REGION#[Pho5.ZRS[1],Pho5.TSS]#

# #X = [len(Pho5._caller_(mol)) for mol in s1]
# X = Molecules.get_region_rval(          s1,*promoter)
# Y = scipy.array([mol._array_ for mol in s1]).transpose()
# correlations = [scipy.stats.pearsonr(X,y)[0] for y in Y[:,]]
# ax.plot(scipy.arange(len(correlations))-Pho5.TSS, Molecules.smooth(correlations,50),color='gold')
# ax.set_ylim([-1,1])

# #X = [len(Pho5._caller_(mol)) for mol in s1]
# X = Molecules.get_region_rval(          s2, *promoter)
# Y = scipy.array([mol._array_ for mol in s2]).transpose()
# correlations = [scipy.stats.pearsonr(X,y)[0] for y in Y[:,]]
# ax.plot(scipy.arange(len(correlations))-Pho5.TSS, Molecules.smooth(correlations,50),color='darkblue')
# ax.set_ylim([-1,1])



##Plotting subsets of the molecule data for different types of analysis
##Below is Rvalues when N-2 is present or not.
# alpha     =0.95
# iterations=1000
# x         =list(range(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS))

# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('R-Value Bootstrapped w/ w/o N-2: {}X'.format(iterations), fontsize=20)

# #subsetting code
# n1     = scipy.array([mol for mol in s6 if int(Pho5._caller_(mol)) in     [0,1,3,5]])
# sans_n1= scipy.array([mol for mol in s6 if int(Pho5._caller_(mol)) not in [0,1,3,5]])

# mew,error,_=Molecules.bootstrap(n1     , Molecules.rvalue, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='gold',linewidth=0.0)
# ax.plot        (x, mew,color='gold',linewidth=2)

# mew,error,_=Molecules.bootstrap(sans_n1, Molecules.rvalue, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='darkblue',linewidth=0.0)
# ax.plot        (x, mew,color='darkblue',linewidth=2)


# n1     = scipy.array([mol for mol in s2 if int(Pho5._caller_(mol)) in     [0,1,3,5]])
# sans_n1= scipy.array([mol for mol in s2 if int(Pho5._caller_(mol)) not in [0,1,3,5]])

# mew,error,_=Molecules.bootstrap(n1     , Molecules.rvalue, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='cyan',linewidth=0.0)
# ax.plot        (x, mew,color='cyan',linewidth=2)

# mew,error,_=Molecules.bootstrap(sans_n1, Molecules.rvalue, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='darkblue',linewidth=0.0)
# ax.plot        (x, mew,color='darkblue',linewidth=2)

##Below is Rvalues for different # of promoter nucleosomes is present or not.
# alpha     =0.95
# iterations=1000
# x         =list(range(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS))

# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('<tata><chd1> Bootstrapped N=0,1,2,3: {}X'.format(iterations), fontsize=20)

# #subsetting code
# subset = s3
# func   = Molecules.rvalue

# n0     = scipy.array([mol for mol in subset if len(Pho5._caller_(mol)) == 0 ])
# n1     = scipy.array([mol for mol in subset if len(Pho5._caller_(mol)) == 1 ])
# n2     = scipy.array([mol for mol in subset if len(Pho5._caller_(mol)) == 2 ])
# n3     = scipy.array([mol for mol in subset if len(Pho5._caller_(mol)) == 3 ])

# # mew,error,_=Molecules.bootstrap(n0     ,func , iterations, alpha)
# # ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='darkblue',linewidth=0.0)
# # ax.plot        (x, mew,color='darkblue',linewidth=2)

# mew,error,_=Molecules.bootstrap(n1     , func, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='lightblue',linewidth=0.0)
# ax.plot        (x, mew,color='darkblue',linewidth=2)

# mew,error,_=Molecules.bootstrap(n2     , func, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='purple',linewidth=0.0)
# ax.plot        (x, mew,color='purple',linewidth=2)

# mew,error,_=Molecules.bootstrap(n3     , func, iterations, alpha)
# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='lightcoral',linewidth=0.0)
# ax.plot        (x, mew,color='darkred',linewidth=2)


##Plotting the rvalue of Bubbles that overlap a specific position
# alpha     =0.95
# iterations=100
# x         =list(range(-Pho5.TSS,Pho5.LENGTH-Pho5.TSS))

# fig,ax = Pho5.build_gene_plot(plt)
# fig.suptitle('Bootestrapped R-Value @ specific position {}X'.format(iterations), fontsize=20)

# uas2= scipy.array(list(Molecules.get_bubbles(s1, *Pho5.UAS2)))
# mew,error,_=Molecules.bootstrap(uas2,
# 					lambda bubbles: Molecules.bubbles_to_molecule(bubbles, 
# 									length=Pho5.LENGTH,
# 									depth =len(uas2))      ,
# 					iterations,
# 					alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='gold',linewidth=0.0)
# ax.plot(x, mew,color='gold',linewidth=2)

# uas2= scipy.array(list(Molecules.get_bubbles(s2, *Pho5.UAS2)))
# mew,error,_=Molecules.bootstrap(uas2,
# 					lambda bubbles: Molecules.bubbles_to_molecule(bubbles, 
# 									length=Pho5.LENGTH,
# 									depth =len(uas2))      ,
# 					iterations,
# 					alpha)

# ax.fill_between(x, error[0,:], error[1,:], alpha=0.30,color='darkblue',linewidth=0.0)
# ax.plot(x, mew,color='darkblue',linewidth=2)





plt.show()