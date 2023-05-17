
import numpy as np
import matplotlib.pyplot as plt

with open("tiempo_DP.dat") as catalog:
    t_DP = []
    for line in catalog:
        column = line.split()
        if not line.startswith('#'): #skipping column labels
            t_DP.append(float(column[0]))

with open("tiempo_RK.dat") as catalog:
    t_RK = []
    for line in catalog:
        column = line.split()
        if not line.startswith('#'): #skipping column labels
            t_RK.append(float(column[0]))
            
with open("energias_RK.dat") as catalog:
    E_RK = []
    for line in catalog:
        column = line.split()
        if not line.startswith('#'): #skipping column labels
            E_RK.append(float(column[0]))

with open("energias_DP.dat") as catalog:
    E_DP = []
    for line in catalog:
        column = line.split()
        if not line.startswith('#'): #skipping column labels
            E_DP.append(float(column[0]))

E0 = E_DP[0]

E_DP = np.array(E_DP)
E_RK = np.array(E_RK)

plt.plot(t_DP, abs(E_DP/E0), label = r'$E_{DP}$')
plt.plot(t_RK[:-1], abs(E_RK/E0), label = r'$E_{RK 4}$')
plt.minorticks_on()
plt.tick_params(axis="both", length=5, direction="in", bottom=True, top=True, left=True, right=True)
plt.tick_params(which='minor', length=2, axis="both", direction="in", bottom=True, top=True, left=True, right=True)
plt.xlabel(r'$t$', fontsize=12)
plt.ylabel(r'$|E/E_0|$', fontsize=12)
plt.axhline(y=E0/E0, lw=2, ls='dashed', alpha=0.5)
plt.legend(loc='lower right')
plt.savefig('energias.png')
plt.show()
