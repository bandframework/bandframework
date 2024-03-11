import matplotlib.pyplot as plt
import numpy as np
import os
from pylab import *
from matplotlib.colors import LogNorm
plt.rc('text', usetex=False)

# you will likely wish to edit these parameters
ParsToPlot=[1,2,3,4,5,0]    #Choose which model parameters to plot, and their order
npanels=len(ParsToPlot)
nbins=30 # Higher gives better resolution but can be noisier
mcmcfilename='trace.txt' # Created by MCMC program (look in mcmc_trace/..)
limitsfilename='modelpar_info.txt' #same as used by Smooth Emulator, but needs extra entry for parname for plot
# This parameter name can be "latex-like", e.g. $\alpha_s$
outputfilename='posterior.pdf'

#you might wish to edit these
ticklabelsize=20.0/(1+npanels**0.25)
labelsize=1.4*ticklabelsize
marginpadding=0.27 # Margin outside plots. Adjust to fit in labels and ticklabels
ThetaMax=1.2 # maximum range of plots

###########################################
#hopefully, everything below will just work
squaresize=1.0
totalxsize=2.0*marginpadding+npanels*squaresize
totalysize=totalxsize
plt.figure(figsize=(totalxsize,totalysize))
fig = plt.figure(1)
marginpadding=marginpadding/totalxsize
squaresize=squaresize/totalxsize
counts1d = np.ndarray((nbins),float)
xarray = np.ndarray((nbins),float)

mcmcdata=np.loadtxt(mcmcfilename,skiprows=0,unpack=True)
npts = mcmcdata[0].size
print('Number of points in trace =', npts)

limitsdata=[[],[]]
iline=0
npars=0
for line in open ('modelpar_info.txt'):
  if len(line) > 3:
    npars=npars+1
    limitsdata.append([])
    tempdata=list(map(str,line.split()))
    for iword in range(0,len(tempdata)):
      limitsdata[iline].append(tempdata[iword])
    iline=iline+1

for ipar in range (0,npars):
  nstrings = len(limitsdata[ipar])
  for istring in range(3,nstrings):
    limitsdata[ipar][2] = str(limitsdata[ipar][2])+' '+str(limitsdata[ipar][istring])
  #print(limitsdata[ipar])

xmin=-1.0
xmax=1.0
ymin=-1.0
ymax=1.0

for ipanel in range (0,npanels):
  for jpanel in range (0,npanels):
    ipar=ParsToPlot[ipanel]
    jpar=ParsToPlot[jpanel]
    
    ax = fig.add_axes([marginpadding+ipanel*squaresize,marginpadding+(npanels-jpanel-1)*squaresize,squaresize,squaresize])
    ax.set_xticks([-1,0,1], minor=False)
    ax.set_yticks([-1,0,1],minor=False)
    ax.tick_params(axis='y',left=True,right=True,labelleft=False,labelright=False,labelsize=ticklabelsize)
    ax.tick_params(axis='x',top=True,bottom=True,labeltop=False,labelbottom=False,labelsize=ticklabelsize)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if ipanel == 0:
      ax.tick_params(axis='y',labelleft=True,labelright=False,labelsize=ticklabelsize)
      ax.set_yticklabels([-1,0,1],minor=False)
      
    if ipanel == npanels-1:
      ax.yaxis.set_label_position('right')
      ax.set_ylabel(limitsdata[jpanel][2],size=labelsize, family='serif', rotation=-90, va='bottom',labelpad=0.1*labelsize)
    
    if jpanel == 0:
      ax.xaxis.set_label_position('top')
      ax.set_xlabel(limitsdata[ipar][2],size=labelsize, family='serif',labelpad=0.2*labelsize)
       
    if jpanel == npanels-1:
      ax.tick_params(axis='x',labelbottom=True,labeltop=False,labelsize=ticklabelsize)
      ax.set_xticklabels([-1,0,1],minor=False)
      
    if ipanel != jpanel:
      counts, xbins, ybins, Image =  ax.hist2d(mcmcdata[ipar],mcmcdata[jpar], bins=nbins, norm=LogNorm(), cmap=plt.cm.YlGnBu)
      counts = counts.transpose()
      maxcounts = np.nanmax(counts)
      level3=maxcounts*exp(-0.5)
      level2=maxcounts*exp(-2.0)
      level1=maxcounts*exp(-4.5)
      ax.contour(counts,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
       linewidths=2.5,levels=(level1,level2,level3),colors=('darkred','darkgreen','darkblue'))
      plt.xlim(-ThetaMax,ThetaMax)
      plt.ylim(-ThetaMax,ThetaMax)
      
    if ipanel==jpanel:
      for ibin in range(0,nbins):
        maxcount=0
        counts1d[ibin]=0.0;
        xarray[ibin] = xmin+(ibin+0.5)*(xmax-xmin)/float(nbins)
        
      for i in range(0,npts):
        x=mcmcdata[ipar][i]
        ibin=int(np.floor(nbins*(x-xmin)/(xmax-xmin)))
        if ibin >= nbins or ibin < 0:
          print ('Outside (-1,1) interval, Theta = ', x)
        else:
          counts1d[ibin] = counts1d[ibin]+1.0
          if counts1d[ibin] > maxcount:
            maxcount = counts1d[ibin]

      counts1d = counts1d/float(maxcount)      
      ax.set_yticks([0,0.5,1.0], minor=False)
      ax.tick_params(axis='y',labelleft=True,labelright=False,labelsize=ticklabelsize)
      plt.ylim(-0.1,1.1)
      plt.xlim(-ThetaMax,ThetaMax)
      plt.plot(xarray,counts1d,linestyle='-',linewidth=3,color='r')
    
#  plt.show()
plt.savefig(outputfilename,format='pdf')
plt.show()
plt.close()
# if you have Mac OS and want to see pdf file, comment out previous two lines and uncomment line below
#os.system("open -a Preview "+outputfilename);

quit()
