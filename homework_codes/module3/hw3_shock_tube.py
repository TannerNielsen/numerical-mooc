#!/usr/bin/python

# MAE6286 Practical Numerical Methods with Python
#
# Module 3
# 
# The program solves Sod's Shock Tube problem (1D Euler Equations)
#
#        |  rho  |  +       |    rho*u    |    |0|
#   d/dt | rho*u |  +  d/dx | rho*u^2 + p |  = |0| 
#        |rho*e_T|  +       |(rho*e_T+p)*u|    |0|
#
#  This code uses the Richtmyer method, a two-step prediction correction method, to advance the solution in time.
#
#  Tanner Nielsen 2-16-2015

import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def computeF(u,nx,gamma):
  flux = np.zeros((nx,3))
  flux[:,0] = u[:,1]
  flux[:,1] = u[:,1]**2/u[:,0] + (gamma-1.0)*(u[:,2]-u[:,1]**2/(2.0*u[:,0]))
  flux[:,2] = (u[:,2] + (gamma-1.0)*(u[:,2]-u[:,1]**2/(2.0*u[:,0])))*(u[:,1]/u[:,0])
  return flux


def richtmyer(u,dx,dt,nx,gamma):
  un = np.copy(u)
  ustar = np.copy(u)
  F = computeF(u,nx,gamma)
  ustar[:-1,:] = 0.5*(un[1:,:]+un[:-1,:]) - dt/(2.0*dx)*(F[1:,:] - F[:-1,:])
  Fstar = computeF(ustar,nx,gamma)
  u[1:,:] = un[1:,:] - dt/dx*(Fstar[1:,:] - Fstar[:-1,:])
  return u


def main():

  L = 20.0
  nx = 81
  dx = L/(nx-1)
  print 'dx: ', dx
  dt = 0.0002
  nt = 51
  gamma = 1.4

  # initialize density
  rho = np.ones(nx)*1.0
  rho[nx/2:] = 0.125
  # initialize velocity
  v = np.zeros(nx)
  # initialize pressure
  p = np.ones(nx)*100000.0
  p[nx/2:] = 10000.0

  # initialize solution vector
  u = np.zeros((nx,3))

  for i in range(1,nt):
    u[:,0] = rho
    u[:,1] = rho*v
    u[:,2] = rho*(p/((gamma-1.0)*rho)+0.5*v**2)
    print 'Time: ', dt*i
    u = richtmyer(u,dx,dt,nx,gamma)
    rho = u[:,0]
    v = u[:,1]/u[:,0]
    p = (gamma-1.0)*(u[:,2] - u[:,1]**2/(2.0*u[:,0]))

    # plot rho
    plt.figure(1)
    plt.plot(np.linspace(-10,10,nx), rho, color='red', ls='-', lw=3)
    plt.title(r'Density',fontsize=18)
    pylab.savefig('./plots/shock_tube/density_'+str(i)+'.png')
    plt.close(1);
    # plot velocity
    plt.figure(2)
    plt.plot(np.linspace(-10,10,nx), v, color='red', ls='-', lw=3)
    plt.title(r'Velocity',fontsize=18)
    pylab.savefig('./plots/shock_tube/velocity_'+str(i)+'.png')
    plt.close(2);
    # plot pressure
    plt.figure(3)
    plt.plot(np.linspace(-10,10,nx), p, color='red', ls='-', lw=3)
    plt.title(r'Pressure',fontsize=18)
    pylab.savefig('./plots/shock_tube/pressure_'+str(i)+'.png')  
    plt.close(3);

  print 'The velocity at x=2.5m is: ', v[int(12.5/dx)]
  print 'The pressure at x=2.5m is: ', p[int(12.5/dx)]
  print 'The density at x=2.5m is: ', rho[int(12.5/dx)]


if __name__ == '__main__':
  main()
