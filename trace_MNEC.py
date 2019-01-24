import matplotlib.pyplot as plt
import numpy as np
from math import *
from pylab import *

flowRusa = np.loadtxt("EF1D_Rusanov.txt")
flowRelax = np.loadtxt("EF1D_Relaxation.txt")

RusaNorme = np.loadtxt("EF1D_Rusanov_norme.txt")

flowRusa = np.transpose(flowRusa)
flowRelax = np.transpose(flowRelax)

RusaNorme = np.transpose(RusaNorme)

fig,ax = plt.subplots()

# plt.plot(RusaNorme[0],RusaNorme[1],
# color='green',linewidth=2, label='Rusanov')

plt.plot(flowRelax[0],flowRelax[2],
color='red', marker='', linestyle='-',linewidth=2, markersize=4, label='rho_eq')

plt.plot(flowRelax[0],flowRelax[1],
color='blue', marker='', linestyle='-',linewidth=2, markersize=4, label='Relaxation - Dirichlet')

plt.ylabel("densite")
plt.xlabel("x(m)")
plt.legend(loc=0)

# plt.axis([-0.06, 0.06, 0.01, 0.7])


filename="rho_relax.png"
fig.savefig(filename)
