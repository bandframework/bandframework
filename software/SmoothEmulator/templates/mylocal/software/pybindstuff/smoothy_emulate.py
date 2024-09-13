import emulator_smooth
import numpy as np
smoothmaster=emulator_smooth.emulator_smooth()

NPars=smoothmaster.GetNPars()
NObs=smoothmaster.GetNObs()

theta=np.zeros(NPars,dtype='float')
Y=np.zeros(NObs,dtype='float')
SigmaY=np.zeros(NObs,dtype='float')

smoothmaster.TuneAllY()
print('tuning completed')

print('theta=',theta)

for iY in range(0,NObs):
  #Y[iY]=a.GetYOnlyPython(iY,theta)
  #print('Y[',iY,']=',Y[iY])
  Y[iY],SigmaY[iY]=smoothmaster.GetYSigmaPython(iY,theta)
  print('Y[',iY,']=',Y[iY],' SigmaY=',SigmaY[iY])
  