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

plt.plot(RusaNorme[0],RusaNorme[1],
color='green',linewidth=2, label='Rusanov')

# plt.plot(flowRelax[0],flowRelax[2],
# color='red', marker='o', linestyle='',linewidth=2, markersize=4, label='rho_eq')


plt.ylabel("erreur en norme L2 de la densite")
plt.xlabel("T(s)")
plt.legend(loc=0)

# plt.axis([-0.06, 0.06, 0.01, 0.7])


filename="Norme_rho_rusanov.png"
fig.savefig(filename)
