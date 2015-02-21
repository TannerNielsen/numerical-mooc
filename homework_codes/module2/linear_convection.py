#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pylab

c = 1.

def du_dx(u,dx):

#  u_x = (u[1:] - u[0:-1])/dx
#  u_x[0] = u[0]

  u_x = np.zeros(len(u))

  for i in range(len(u)):
    # backward differencing at right boundary
    if i == len(u)-1:
      u_x[i] = (u[i]-u[i-1])/dx
    # forward differencing everywhere else
    else:
      u_x[i] = (u[i+1]-u[i])/dx

  return u_x

def f(u,dx):

 return -c*du_dx(u,dx) 

def euler_step(u,t,dt):

  return u + dt*f(u,t)

def main():

  nx = 5
  dx = 2./(nx-1)
  nt = 25
  dt = 0.05
#  c = 1.0

  u = np.ones(nx)
  u[.5/dx:1/dx+1] = 2
#  print .5/dx,1/dx
#
#  print u

#  u = np.zeros(nx)
#  v = np.zeros(nx)

  x = np.linspace(0,2,nx)
#  t = np.linspace(0,dt*nt,dt)

#  print t

  for n in range(nt):
    u[n+1] = euler_step(u[n],n,dt)

#  u = 2.*x[:]
#  v = du_dx(u,dx)

#  print len(x),len(u),len(v)

  plt.plot(np.linspace(0,2,nx),u)
#  plt.plot(x,v)
  plt.ylim(0,2.5);
  plt.show()


if __name__ == '__main__':
  main()
