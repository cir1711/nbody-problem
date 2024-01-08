import matplotlib.pyplot as plt
import numpy as np

def plotea(file_name, png_name):
  with open(file_name) as f:
    l  = [i.split("\t") for i in f.read().split("\n")]
    l2 = [[float(j) for j in i if j != ""] for i in l]

  X = np.zeros([int(len(l2[0])/2), len(l2)-1])
  Y = np.zeros([int(len(l2[0])/2), len(l2)-1])


  for i in range(int(len(l2[0])/2)):
    for j in range(len(l2)-1):
      X[i][j] = l2[j][2*i]
      Y[i][j] = l2[j][2*i+1]
    

  plt.figure()
  plt.minorticks_on()
  plt.tick_params(axis="both", length=5, direction="in", bottom=True, top=True, left=True, right=True)
  plt.tick_params(which='minor', length=2, axis="both", direction="in", bottom=True, top=True, left=True, right=True)
  plt.xlabel(r'$x$', fontsize=12)
  plt.ylabel(r'$y$', fontsize=12)
  for i in range(len(X)):
    plt.plot(X[i],Y[i])
  plt.savefig(png_name)
  plt.show()


plotea("./../results/posiciones_DP.dat", "./../plots/orbita_DP")
plotea("./../results/posiciones_RK.dat", "./../plots/orbita_RK")