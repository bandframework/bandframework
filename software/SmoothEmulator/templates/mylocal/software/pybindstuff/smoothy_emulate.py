import numpy as np
import sys
sys.path.insert(0,"/Users/CarlosSmith/mylocal/software/pybindstuff")
import emulator_smooth
smoothmaster=emulator_smooth.emulator_smooth()

smoothmaster.TuneAllY()
print('tuning completed')

NPars=smoothmaster.GetNPars()
NObs=smoothmaster.GetNObs()

theta=np.zeros(NPars,dtype='float')
X=np.zeros(NPars,dtype='float')

X[0]=210.0
X[1]=0.1715
X[2]=0.75
X[3]=0.55
X[4]=1.4
X[5]=22.5

print('X=',X)
theta=smoothmaster.GetThetaFromX(X)
print('theta=',theta)

Y=np.zeros(NObs,dtype='float')
SigmaY=np.zeros(NObs,dtype='float')

for iY in range(0,NObs):
  #Y[iY]=a.GetYOnlyPython(iY,theta)
  #print('Y[',iY,']=',Y[iY])
  Y[iY],SigmaY[iY]=smoothmaster.GetYSigmaPython(iY,theta)
  print('Y[',iY,']=',Y[iY],' SigmaY=',SigmaY[iY])
  
