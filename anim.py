import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc


def anim(file_name, time_file_name, video_name):
  with open(file_name) as f:
    l  = [i.split("\t") for i in f.read().split("\n")]
    l2 = [[float(j) for j in i if j != ""] for i in l]

  X = np.zeros([int(len(l2[0])/2), len(l2)-1])
  Y = np.zeros([int(len(l2[0])/2), len(l2)-1])
  x = [[] for i in range(len(X))]
  y = [[] for i in range(len(X))]


  for i in range(int(len(l2[0])/2)):
    for j in range(len(l2)-1):
      X[i][j] = l2[j][2*i]
      Y[i][j] = l2[j][2*i+1]

  with open(time_file_name) as catalog:
    t = []
    for line in catalog:
        column = line.split()
        if not line.startswith('#'): #skipping column labels
            t.append(float(column[0]))


  n = 0
  while n*0.01<=max(t):
    for i in range(len(t)):
      if t[i] >= n*0.01:
        for j in range(len(X)):
          x[j].append(X[j][i])
          y[j].append(Y[j][i])
        n += 1



  M = len(x[0])

  fig1 = plt.figure()

  l1, = plt.plot([], [], 'o--',lw=2, markersize=8, markevery=[-1])
  l2, = plt.plot([], [], 'o--',lw=2, markersize=8, markevery=[-1])
  l3, = plt.plot([], [], 'o--',lw=2, markersize=8, markevery=[-1])
  def update_line(num):
      num +=1
      l1.set_data(x[0][:num],y[0][:num])
      l2.set_data(x[1][:num],y[1][:num])
      l3.set_data(x[2][:num],y[2][:num])
      return l1, l2, l3,
  ymin = min([min(y[0]),min(y[1]),min(y[2])])
  ymax = max([max(y[0]),max(y[1]),max(y[2])])
  xmin = min([min(x[0]),min(x[1]),min(x[2])])
  xmax = max([max(x[0]),max(x[1]),max(x[2])])
  plt.xlim(xmin-0.1,xmax+0.1)
  plt.ylim(ymin-0.1,ymax+0.1)
  plt.minorticks_on()
  plt.tick_params(axis="both", length=5, direction="in", bottom=True, top=True, left=True, right=True)
  plt.tick_params(which='minor', length=2, axis="both", direction="in", bottom=True, top=True, left=True, right=True)
  plt.xlabel(r'$x$', fontsize=12)
  plt.ylabel(r'$y$', fontsize=12)
  line_ani = animation.FuncAnimation(fig1, update_line, frames = M,
      interval=10, blit=True) #el interval ves cambiando como te guste, pero con 10 ya va bien
  line_ani.save(video_name)
  

anim("posiciones_DP.dat", "tiempo_DP.dat","anim_DP.mp4")