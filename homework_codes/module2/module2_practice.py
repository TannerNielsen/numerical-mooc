#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def main():

  nx = 41
  dx = 2./(nx-1)
  print 'dx: ', dx
  nt = 25
  dt = 0.02
  c = 1.0

  u = np.ones(nx)      #numpy function ones()
  u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s
  print(u)

  plt.plot(np.linspace(0,2,nx), u, color='red', ls='--', lw=3)
  plt.ylim(0,2.5);

  for n in range(1,nt):  
    un = u.copy() 
    for i in range(1,nx): 
        u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1])  # nonlinear wave equation: du/dt + u*du/dx = 0
        #u[i] = un[i]-c*dt/dx*(un[i]-un[i-1])  # linear wave equation: du/dt + c*du/dx = 0

  
  plt.plot(np.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
  plt.ylim(0,2.5);
  pylab.show()

  


if __name__ == '__main__':
  main()
