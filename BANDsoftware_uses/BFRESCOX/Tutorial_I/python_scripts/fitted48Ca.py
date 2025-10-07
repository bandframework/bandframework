import os
import matplotlib.pyplot as plt
import numpy as np
import pylab

#os.chdir("/Users/ozgesurer/Desktop/Frescox/experiment/")

# Execute frescox
os.system("frescox < frescox_inputs/48Ca.in frescox_inputs/search-el.in frescox_inputs/sfresco.in > frescox_outputs/48Ca.out")

# Read outputs
files = ['frescox_outputs/48Ca.out']
outs = []
for file in files:
    with open(file) as f:
        content = f.readlines()
        outs.append(content)
        
slist = []
for out in outs:
    sigma_omega_ratio = [] 
    for idline, line in enumerate(out):
        if ('X-S' in line):
            sigma_omega_ratio.append(float(line.split()[4]))

    slist.append(sigma_omega_ratio)

# Observed data
data = np.array([[25.5, 1243], 
                 [30.6, 887.7],               
                 [40.8, 355.5],               
                 [50.9, 111.5],               
                 [61, 26.5],             
                 [71.1, 10.4],           
                 [76.2, 8.3],             
                 [81.2, 7.3],             
                 [91.2, 17.2],             
                 [101.2, 37.6],              
                 [111.1, 48.7],                
                 [121, 38.9],              
                 [130.9, 32.4],               
                 [140.8, 36.4],               
                 [150.6, 61.9]])



fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1)

line, = ax.plot(slist[0], color='red', label='Best Fit Parametrization')
ax.scatter(data[:,0], data[:,1], label='Elastic Scattering Data')
plt.xlabel(r'$\theta$ (deg)')
plt.ylabel(r'd$\sigma/d\Omega$')
#plt.legend(bbox_to_anchor=(0.5, -0.2), loc='center', borderaxespad=0.)
plt.legend()
ax.set_yscale('log')

pylab.show()

