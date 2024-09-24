import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from pylab import *

# User will likely wish to adjust many of these
writetitle=True
topmargin=0.01   # increase if you want title
if writetitle:
  topmargin=0.07
bottommargin=0.15 # adjust to fit x-axis labels
leftmargin=0.36    # adjust to fit y-axis labels
thetamax=3.5      # sets max display range of y axes
myfontsize=9      # base of all fonts
barwidth=0.5      # larger number stretches figure horizontally
panelheight=0.75  # larger values stretches figure vertically

# User can hopefully leave all this alone
Ynames = []
NObs=0
iline=0
for line in open('../observable_info.txt'):
  if len(line) > 1:
    NObs+=1
    tempdata=list(map(str,line.split()))
    Ynames.append(tempdata[1])
    for iword in range(2,len(tempdata)):
      Ynames[iline]=Ynames[iline]+' '+tempdata[iword]
    iline+=1

Xnames=[]
NPars=0
iline=0
for line in open ('../modelpar_info.txt'):
  if len(line) > 2:
    tempdata=list(map(str,line.split()))
    NPars=NPars+1
    Xnames.append(tempdata[2])
    for iword in range(3,len(tempdata)):
      Xnames[iline]=Xnames[iline]+' '+tempdata[iword]
    iline+=1
    
RPdata = np.loadtxt("ResolvingPower.txt",skiprows=0,unpack=True)
Yi = arange(0,NObs-0.01,1)

plt.figure(figsize=(leftmargin+NPars*barwidth,bottommargin+NObs*panelheight))
fig = plt.figure(1)

plotwidth=(1.0-leftmargin)-0.01
plotheight=1.0-bottommargin-topmargin

Npanels=NPars
for ipanel in range (0,Npanels):
  ax = fig.add_axes([leftmargin,bottommargin+ipanel*plotheight/Npanels,
  plotwidth,plotheight/Npanels])
  ax.set_yticks(np.arange(-6,6,2), minor=False)
  ax.set_yticks(np.arange(-6,6,1), minor=True)
  ax.set_yticklabels(np.arange(-6,6,2), minor=False, family='serif', size=myfontsize)

  xi = RPdata[ipanel]
  xxi = 0.5*(np.abs(xi)+xi)
  xxxi= 0.5*(-np.abs(xi)+xi)
    
  plt.ylabel(Xnames[ipanel],fontsize=1.3*myfontsize,position=(0.1,0.35),rotation=0.0,ha='right')
  plt.bar(Yi,xxi,0.8,color='r')
  plt.bar(Yi,xxxi,0.8,color='b')
  
  ax.set_xticks(np.arange(0,NObs-0.1,1),minor=False)
  if ipanel == 0:
    ax.set_xticklabels(Ynames, minor=False, size=1.3*myfontsize,rotation=90)
  else:
    ax.set_xticklabels([], minor=False)
  
  if writetitle:
    ax.set_title('$\partial \langle\langle \\theta_i\\rangle\\rangle/\partial Y^{\\rm exp}_a \langle \delta Y_a^2\\rangle^{1/2}$',fontsize=1.3*myfontsize)

  plt.xlim(-0.5,NObs-0.5)
  plt.ylim(-thetamax,thetamax)

plt.savefig('RP.pdf',format='pdf')

plt.show()
plt.close()
# if you have Mac OS and want to see pdf file, comment out previous two lines and uncomment line below
#os.system("open -a Preview RP.pdf") # syntax works for Mac OS

quit()
