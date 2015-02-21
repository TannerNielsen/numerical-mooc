#!/usr/bin/python

# MAE6286 Practical Numerical Methods with Python
# Module 4:
#
# This code explicitly solves 2D reaction-diffusion equations (Gray-Scott Model)
#
#  du/dt = D_u*(d^2u/dx^2 + d^2u/dy^2) - u*v^2 + F*(1-u)
#  dv/dt = D_v*(d^2v/dx^2 + d^2v/dy^2) + u*v^2 - (F+k)*v
#
#  where u and v are two generic chemical species.
#
#  Tanner Nielsen 2-20-2015

import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys
import matplotlib.cm as cm
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def main():

  # input plot output name and constants
  #plotname = 'bacteria_1_'
  #plotname = 'bacteria_2_'
  #plotname = 'coral_'
  #plotname = 'fingerprint_'
  #plotname = 'spirals_'
  #plotname = 'spirals_dense_'
  #plotname = 'spirals_fast_'
  #plotname = 'unstable_'
  #plotname = 'worms_1_'
  #plotname = 'worms_2_'
  plotname = 'zebrafish_'

  #Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1
  #Du, Dv, F, k = 0.00014, 0.00006, 0.035, 0.065 # Bacteria 2
  #Du, Dv, F, k = 0.00016, 0.00008, 0.060, 0.062 # Coral
  #Du, Dv, F, k = 0.00019, 0.00005, 0.060, 0.062 # Fingerprint
  #Du, Dv, F, k = 0.00010, 0.00010, 0.018, 0.050 # Spirals
  #Du, Dv, F, k = 0.00012, 0.00008, 0.020, 0.050 # Spirals Dense
  #Du, Dv, F, k = 0.00010, 0.00016, 0.020, 0.050 # Spirals Fast
  #Du, Dv, F, k = 0.00016, 0.00008, 0.020, 0.055 # Unstable
  #Du, Dv, F, k = 0.00016, 0.00008, 0.050, 0.065 # Worms 1
  #Du, Dv, F, k = 0.00016, 0.00008, 0.054, 0.063 # Worms 2
  Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.060 # Zebrafish

  # initialization
  uvinitial = np.load('./uvinitial.npz')
  U = uvinitial['U']
  V = uvinitial['V']

  # plot initialization
  plt_count = 1  # initialize plot counter
  fig = plt.figure(figsize=(8,5))
  plt.suptitle('Gray-Scott Model', fontsize=18)
  plt.subplot(121)
  plt.imshow(U, cmap = cm.RdBu)
  plt.xticks([]), plt.yticks([]);
  plt.title('Species U')
  plt.subplot(122)
  plt.imshow(V, cmap = cm.RdBu)
  plt.xticks([]), plt.yticks([]);
  plt.title('Species V')
  pylab.savefig('./plots/'+plotname+str(plt_count)+'.png')
  plt.close()

  nx = 192
  ny = 192
  print 'Du, Dv, F, k = ', Du, Dv, F, k
  T = 8000
  Lx = 5.0
  Ly = 5.0
  dx = Lx/(nx-1)
  dy = Ly/(ny-1)
  dt = .9 * min(dx,dy)**2 / (4*max(Du,Dv))
  nt = int(T/dt)
  print 'dt = ', dt
  print 'Number of time steps = ', nt
  print 'dx, dy = ', dx, dy
  print 'nx, ny = ', nx, ny
  print 'Domain size (x,y) = ', Lx, Ly

  # plot frequency and counter
  plot_freq = 100
  pcount = 0

  # initialize u and v matrices
  u = U.copy()
  v = V.copy()

  for n in range(nt):
    print 'Time step, time = ', n+1, (n+1)*dt
    pcount = pcount + 1

    un = u.copy()
    vn = v.copy()

    # explicitly solve u equation
    u[1:-1,1:-1] = un[1:-1,1:-1] + dt*(Du*((un[1:-1,2:] - 2.0*un[1:-1,1:-1] + un[1:-1,:-2])/dx**2 +
                                           (un[2:,1:-1] - 2.0*un[1:-1,1:-1] + un[:-2,1:-1])/dy**2) -
                                       un[1:-1,1:-1]*vn[1:-1,1:-1]**2 + F*(1.0-un[1:-1,1:-1]))
    # enforce Neumann BCs for u equation
    u[-1,:] = u[-2,:]  # top
    u[:,-1] = u[:,-2]  # right
    u[0,:] = u[1,:]    # bottom
    u[:,0] = u[:,1]    # left

    # explicitly solve v equation
    v[1:-1,1:-1] = vn[1:-1,1:-1] + dt*(Dv*((vn[1:-1,2:] - 2.0*vn[1:-1,1:-1] + vn[1:-1,:-2])/dx**2 +
                                           (vn[2:,1:-1] - 2.0*vn[1:-1,1:-1] + vn[:-2,1:-1])/dy**2) +
                                       un[1:-1,1:-1]*vn[1:-1,1:-1]**2 - (F+k)*vn[1:-1,1:-1])

    # enforce Neumann BCs for v equation
    v[-1,:] = v[-2,:]  # top
    v[:,-1] = v[:,-2]  # right
    v[0,:] = v[1,:]    # bottom
    v[:,0] = v[:,1]    # left

    # plot results
    if pcount == plot_freq:
      print 'Plotting iteration: ', n+1
      plt_count = plt_count + 1
      fig = plt.figure(figsize=(8,5))
      plt.suptitle('Gray-Scott Model', fontsize=18)
      plt.subplot(121)
      plt.imshow(u, cmap = cm.RdBu)
      plt.xticks([]), plt.yticks([]);
      plt.title('Species U')
      plt.subplot(122)
      plt.imshow(v, cmap = cm.RdBu)
      plt.xticks([]), plt.yticks([]);
      plt.title('Species V')
      pylab.savefig('./plots/'+plotname+str(plt_count)+'.png')
      plt.close()
      pcount = 0
    
  print 'The answer is: ', u[100,::40]


if __name__ == '__main__':
  main()
