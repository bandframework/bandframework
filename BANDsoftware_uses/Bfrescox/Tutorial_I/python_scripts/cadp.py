import os
import matplotlib.pyplot as plt

#os.chdir("/Users/ozgesurer/Desktop/Frescox/experiment/")
#os.system("../source/frescox < D1048cadp.in > D1048cadp.out")
#os.system("../source/frescox < LH1048cadp.in > LH1048cadp.out")

os.system("frescox < frescox_inputs/D1048cadp.in > frescox_outputs/D1048cadp.out")
os.system("frescox < frescox_inputs/LH1048cadp.in > frescox_outputs/LH1048cadp.out")

# Read the output files
files = ['frescox_outputs/LH1048cadp.out', 'frescox_outputs/D1048cadp.out']
outs = []
for file in files:
    with open(file) as f:
        content = f.readlines()
        outs.append(content)

# Get the corresponding outputs
Rlist, slist = [], []
for out in outs:
    Rutherford = [] 
    sigma_omega_ratio = [] 
    for idline, line in enumerate(out):
        if '/R' in line:
            Rutherford.append(float(line.split()[3]))
        if (idline > 575) and ('X-S' in line):
            sigma_omega_ratio.append(float(line.split()[4]))
    
    Rlist.append(Rutherford)
    slist.append(sigma_omega_ratio)

# Plot the outputs
i = 0
fig, axes = plt.subplots(1, 2, figsize=(9,4))
for values in [Rlist, slist]:
    axes[i].plot(values[0][0:80], linestyle='dotted', c='blue', label='LH (global)')
    axes[i].plot(values[1][0:80], linestyle='dashed', c='green', label='D (global)')
    i += 1

#plt.legend(bbox_to_anchor=(-0.5, -0.25), loc='center', borderaxespad=0)
fig.subplots_adjust(bottom=0.35)
plt.legend(bbox_to_anchor=(0, -0.4), loc='center')
axes[0].set_xlabel(r'$\theta$ (deg)')
axes[0].set_ylabel(r'd$\sigma/d\Omega$ (ratio to Rutherford)')
axes[1].set_xlabel(r'$\theta$ (deg)')
axes[1].set_ylabel(r'd$\sigma/d\Omega$ (mb/srad)')
plt.show()