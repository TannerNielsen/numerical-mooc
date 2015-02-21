#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pylab

def u_initial(nx):
  u = np.ones(nx)
  u[nx/2:] = 0.0
  return u


def computeF(u):
  return u*u/2.0


def maccormack(u, dt, dx):
  un = u.copy()
  ustar = u.copy()    
  F = computeF(u)
  ustar[:-1] = un[:-1] - dt/dx*(F[1:] - F[:-1])
  Fstar = computeF(ustar)
  u[1:] = 0.5*(un[1:] + ustar[1:] - dt/dx*(Fstar[1:] - Fstar[:-1]))
  
  return u

def maccormack_damping(u, dt, dx, nx, epsilon):
  un = u.copy()
  ustar = u.copy()
  F = computeF(u)
  ustar[:-1] = un[:-1] - dt/dx*(F[1:] - F[:-1])
  ustar[1:nx-2] = ustar[1:nx-2] + epsilon*(un[2:nx-1] - 2.0*un[1:nx-2] + un[0:nx-3])
  Fstar = computeF(ustar)
  u[1:] = 0.5*(un[1:] + ustar[1:] - dt/dx*(Fstar[1:] - Fstar[:-1]))

  return u


def main():


  nx = 81
  nt = 70
  dx = 4.0/(nx-1)
  epsilon = 0.1

  u = u_initial(nx)
  sigma = 0.5
  dt = sigma*dx

  for i in range(1,nt):
    print 'Iter: ', i
    #u = maccormack(u,dt,dx)
    u = maccormack_damping(u, dt, dx, nx, epsilon)
    plt.plot(np.linspace(0,4,nx), u, color='red', ls='-', lw=3)
    plt.title(r"Burger's Equation", fontsize=18)
    plt.ylabel(r'Velocity (m/s)', fontsize=18)
    plt.xlabel(r'X (m)', fontsize=18)
    plt.ylim(0,1.2);
    #pylab.savefig('./plots/burgers_maccormack_'+str(i)+'.png')
    pylab.savefig('./plots/damping/burgers_maccormack_'+str(epsilon)+'_'+str(i)+'.png')
    plt.close();



if __name__ == '__main__':
  main()
